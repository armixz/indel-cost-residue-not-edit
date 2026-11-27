"""Gate 4: is insertion tolerance a property of a flexible SEGMENT rather than of one residue?
See ../RESULTS.md.

run_len   = length of the contiguous run of positions with pLDDT below the protein's own median
            that contains this position (0 if the position itself is above it)
run_slack = min(residues of that run to the left, residues to the right of the position)

Usage: python -I 08_gate4_segments.py <insertion-lab root> <out.json>
"""
import json
import re
import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score
from sklearn.mixture import GaussianMixture

PDB3 = dict(zip(["ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", "ILE", "LEU", "LYS",
                 "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"], list("ARNDCEQGHILKMFPSTWYV")))
AA1 = dict(zip(["Ala", "Arg", "Asn", "Asp", "Cys", "Glu", "Gln", "Gly", "His", "Ile", "Leu", "Lys",
                "Met", "Phe", "Pro", "Ser", "Thr", "Trp", "Tyr", "Val"], list("ARNDCEQGHILKMFPSTWYV")))


def plddt_profile(pdb):
    out = {}
    for line in open(pdb):
        if line.startswith("ATOM") and line[12:16].strip() == "CA" and line[17:20].strip() in PDB3:
            out[int(line[22:26])] = float(line[60:66])
    return out


def runs(profile, q=0.5):
    """For each residue: run_len and run_slack of the contiguous below-threshold run."""
    idx = sorted(profile)
    thr = float(np.quantile([profile[i] for i in idx], q))
    low = {i: profile[i] < thr for i in idx}
    run_len, run_slack = {}, {}
    i = 0
    while i < len(idx):
        if not low[idx[i]]:
            run_len[idx[i]] = 0; run_slack[idx[i]] = 0; i += 1; continue
        j = i
        while j + 1 < len(idx) and low[idx[j + 1]] and idx[j + 1] == idx[j] + 1:
            j += 1
        L = j - i + 1
        for k in range(i, j + 1):
            run_len[idx[k]] = L
            run_slack[idx[k]] = min(k - i, j - k)
        i = j + 1
    return run_len, run_slack, thr


def cutoff(x):
    g = GaussianMixture(2, random_state=0, n_init=5).fit(np.asarray(x).reshape(-1, 1))
    return float(np.mean(sorted(g.means_.ravel())))


def auroc(lab, v):
    m = np.isfinite(v)
    if m.sum() < 20 or len(np.unique(lab[m])) < 2:
        return None
    return float(roc_auc_score(lab[m], v[m]))


def mavedb_insertions(path, length_filter):
    rows = []
    for r in pd.read_csv(path).to_dict("records"):
        h = str(r.get("hgvs_pro") or ""); s = r.get("score")
        m = re.match(r"p\.([A-Za-z]{3})(\d+)_[A-Za-z]{3}\d+ins([A-Za-z]+)$", h)
        if not m or s is None or not np.isfinite(s):
            continue
        n = len(re.findall(r"[A-Z][a-z]{2}", m.group(3)))
        if n != length_filter:
            continue
        rows.append({"pos": int(m.group(2)), "score": float(s)})
    return pd.DataFrame(rows).groupby("pos")["score"].median().reset_index()


def kir_map(p):
    return p if p <= 115 else (None if p <= 123 else p - 8)


def di_screen_table(root):
    """Domain insertions: positives per the Gate 1 label, pooled over inserted domains."""
    f = pd.read_csv(f"{root}/results/features.csv")
    out = {}
    for prot, g in f.groupby("protein"):
        piv = g.pivot_table(index="res0", columns="domain", values="score")
        keep = piv.notna().any(axis=1)
        out[prot] = pd.DataFrame({"pos": piv.index[keep] + 1,                 # 1-based for the PDB
                                  "y": (piv[keep] > 0).any(axis=1).astype(int).values})
    return out


def main(root, out_path):
    md = f"{root}/data/mavedb"
    ds = f"{root}/data/DI_screen/analysis/input_data"
    rep = {"datasets": {}, "robustness": {}}

    def block(name, pdb, table, label_from, mapper=None, q=0.5):
        prof = plddt_profile(pdb)
        rl, rs, thr = runs(prof, q)
        pos = table["pos"].values
        st = np.array([mapper(int(p)) if mapper else int(p) for p in pos], dtype=object)
        ok = np.array([s is not None and s in prof for s in st])
        t = table[ok].copy()
        st = st[ok]
        t["plddt"] = [prof[s] for s in st]
        t["run_len"] = [rl[s] for s in st]
        t["run_slack"] = [rs[s] for s in st]
        if label_from == "score":
            cut = cutoff(t["score"].values)
            lab = (t["score"].values >= cut).astype(int)
        else:
            lab = t["y"].values.astype(int)
        res = {"n": int(len(t)), "tolerated": int(lab.sum()), "plddt_threshold": round(thr, 2),
               "auroc": {"plddt": auroc(lab, -t["plddt"].values),
                         "run_len": auroc(lab, t["run_len"].values.astype(float)),
                         "run_slack": auroc(lab, t["run_slack"].values.astype(float))}}
        a = res["auroc"]
        res["advantage_run_len_minus_plddt"] = (round(a["run_len"] - a["plddt"], 3)
                                                if a["run_len"] and a["plddt"] else None)
        return res, t

    tables = {}
    # short insertions
    kir = f"{md}/kir21_scores.csv"
    for n in (1, 2, 3):
        tb = mavedb_insertions(kir, n)
        rep["datasets"][f"Kir2.1_ins{n}"], tables[f"Kir2.1_ins{n}"] = block(
            f"Kir2.1_ins{n}", f"{md}/af_P35561.pdb", tb, "score", kir_map)
    tb = mavedb_insertions(f"{md}/scores_00000100-a-1.csv", 1)
    rep["datasets"]["TEM1_ins1"], tables["TEM1_ins1"] = block("TEM1", f"{md}/af_P62593.pdb", tb, "score")
    dh = pd.concat([mavedb_insertions(f"{md}/scores_00001219-a-{i}.csv", 1) for i in (1, 2, 3, 4)])
    dh = dh.groupby("pos")["score"].mean().reset_index()
    rep["datasets"]["DHFR_ins1"], tables["DHFR_ins1"] = block("DHFR", f"{md}/af_P00374.pdb", dh, "score")
    # domain insertions
    for prot, tb in di_screen_table(root).items():
        rep["datasets"][f"{prot}_domain"], tables[f"{prot}_domain"] = block(
            prot, f"{ds}/{prot}_AF.pdb", tb, "y")

    # P1
    short = ["Kir2.1_ins3", "TEM1_ins1", "DHFR_ins1"]
    wins = [k for k in short if (rep["datasets"][k]["advantage_run_len_minus_plddt"] or -1) > 0]
    rep["P1"] = {"short_datasets": short, "run_len_beats_plddt_in": wins, "count": len(wins),
                 "note": "the pre-registered mark said '3 of the 4 short-insertion datasets' but the "
                         "stated counting rule yields only 3 datasets; both readings are reported"}
    # P2
    adv = {k: rep["datasets"][k]["advantage_run_len_minus_plddt"] for k in rep["datasets"]}
    dom = [adv[k] for k in adv if k.endswith("_domain") and adv[k] is not None]
    single = [adv[k] for k in ("Kir2.1_ins1", "TEM1_ins1", "DHFR_ins1") if adv[k] is not None]
    rep["P2"] = {"advantage_by_dataset": adv,
                 "kir_grows_with_size": bool(adv.get("Kir2.1_ins3") is not None and adv.get("Kir2.1_ins1") is not None
                                             and adv["Kir2.1_ins3"] > adv["Kir2.1_ins1"]),
                 "domain_mean_advantage": round(float(np.mean(dom)), 3) if dom else None,
                 "single_max_advantage": round(float(max(single)), 3) if single else None}
    rep["P2"]["domain_above_all_single"] = bool(
        rep["P2"]["domain_mean_advantage"] is not None and rep["P2"]["single_max_advantage"] is not None
        and rep["P2"]["domain_mean_advantage"] > rep["P2"]["single_max_advantage"])
    rep["P2"]["passed"] = bool(rep["P2"]["kir_grows_with_size"] and rep["P2"]["domain_above_all_single"])

    # P3: size sensitivity against slack, within Kir2.1
    a, b = tables["Kir2.1_ins1"].set_index("pos"), tables["Kir2.1_ins3"].set_index("pos")
    common = a.index.intersection(b.index)
    sens = (a.loc[common, "score"] - b.loc[common, "score"]).values
    slack = a.loc[common, "run_slack"].values.astype(float)
    r = float(spearmanr(slack, sens).statistic)
    rng = np.random.default_rng(0)
    bs = [spearmanr(slack[i], sens[i]).statistic for i in (rng.integers(0, len(sens), len(sens)) for _ in range(2000))]
    ci = [round(float(q), 3) for q in np.nanquantile(bs, [0.025, 0.975])]
    rep["P3"] = {"n": int(len(common)), "spearman_slack_vs_size_sensitivity": round(r, 3), "ci95": ci,
                 "passed": bool(r < -0.2 and ci[1] < 0)}

    # robustness: other pLDDT thresholds
    for q in (0.4, 0.6):
        rep["robustness"][f"q{q}"] = {}
        for name, pdb, tb, lab, mp in (("Kir2.1_ins3", f"{md}/af_P35561.pdb", mavedb_insertions(kir, 3), "score", kir_map),
                                       ("TEM1_ins1", f"{md}/af_P62593.pdb", mavedb_insertions(f"{md}/scores_00000100-a-1.csv", 1), "score", None),
                                       ("DHFR_ins1", f"{md}/af_P00374.pdb", dh, "score", None)):
            r2, _ = block(name, pdb, tb, lab, mp, q=q)
            rep["robustness"][f"q{q}"][name] = {"adv": r2["advantage_run_len_minus_plddt"], "auroc": r2["auroc"]}

    rep["verdict"] = {"P1_count_of_3": len(wins), "P2_passed": rep["P2"]["passed"], "P3_passed": rep["P3"]["passed"],
                      "H4_supported": bool(rep["P2"]["passed"])}
    json.dump(rep, open(out_path, "w"), indent=1)
    print(json.dumps(rep, indent=1))


if __name__ == "__main__":
    main(*sys.argv[1:3])
