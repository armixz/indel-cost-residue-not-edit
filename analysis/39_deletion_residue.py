"""Marks D-1 to D-5 of PREREG_deletion_residue.md.

Does the structural pattern of deletion cost depend on which residue is deleted?
All 20 residues are present on the deletion side, unlike insertions (Ala and Gly only).

Usage: indel_env python 39_deletion_residue.py <results_dir> <lib3_csv> <afdir> <out_json>
"""
import csv
import json
import os
import re
import sys

import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
exec(open(os.path.join(HERE, "34_scale.py")).read().split("def designed")[0])  # noqa: S102

# Pace & Scholtz 1998 helix propensity, kcal/mol; higher = worse helix former
PACE = {"A": 0.0, "L": 0.21, "R": 0.21, "M": 0.24, "K": 0.26, "Q": 0.39, "E": 0.40, "I": 0.41,
        "W": 0.49, "S": 0.50, "Y": 0.53, "F": 0.54, "H": 0.61, "V": 0.61, "N": 0.65, "T": 0.66,
        "C": 0.68, "D": 0.69, "G": 1.00, "P": 3.16}
VOL = {"G": 60.1, "A": 88.6, "S": 89.0, "C": 108.5, "D": 111.1, "P": 112.7, "N": 114.1,
       "T": 116.1, "E": 138.4, "V": 140.0, "Q": 143.8, "H": 153.2, "M": 162.9, "I": 166.7,
       "L": 166.7, "K": 168.6, "R": 173.4, "F": 189.9, "Y": 193.6, "W": 227.8}
PAT = re.compile(r"^(.*?)_(ins|del)([A-Z])(\d+)$")


def contrast_by_residue(df, a="H", b="L", nmin=2, col="del"):
    """Per residue: within-domain mean(cost in class a) - mean(cost in class b)."""
    out = {}
    for aa, g in df.groupby("aa"):
        v = [float(x[x["ss"] == a][col].mean() - x[x["ss"] == b][col].mean())
             for _, x in g.groupby("domain")
             if (x["ss"] == a).sum() >= nmin and (x["ss"] == b).sum() >= nmin]
        if len(v) >= 10:
            out[str(aa)] = {"n_domains": len(v), "mean": float(np.mean(v)),
                            "ci95": list(np.quantile(
                                np.asarray(v)[np.random.default_rng(0).integers(
                                    0, len(v), size=(5000, len(v)))].mean(axis=1), [.025, .975]))}
    return out


def build_designed(lib3, afdir):
    csv.field_size_limit(10 ** 7)
    import importlib.util
    sp = importlib.util.spec_from_file_location("d17", os.path.join(HERE, "17_dssp.py"))
    D = importlib.util.module_from_spec(sp)
    sp.loader.exec_module(D)
    wt, rows = {}, []
    for r in csv.DictReader(open(lib3, newline="")):
        nm = r["name"]
        if not nm.startswith("HEEH"):
            continue
        if nm.split("_")[-1].startswith("wt"):
            try:
                wt.setdefault(nm.rsplit("_", 1)[0], []).append(float(r["deltaG"]))
            except ValueError:
                pass
            continue
        m = PAT.match(nm)
        if not m or m.group(2) != "del":
            continue
        try:
            rows.append({"domain": m.group(1), "aa": m.group(3), "pos": int(m.group(4)),
                         "dG": float(r["deltaG"]), "ci": float(r["deltaG_95CI"])})
        except (ValueError, KeyError):
            continue
    v = pd.DataFrame(rows)
    wtm = {k: float(np.mean(x)) for k, x in wt.items()}
    v = v[v["ci"] <= 1.0].copy()
    v["wt"] = v["domain"].map(wtm)
    v = v[v["wt"].notna()]
    v["del"] = v["wt"] - v["dG"]
    ss = []
    for dom in sorted(v["domain"].unique()):
        p = os.path.join(afdir, dom + ".pdb")
        if not os.path.exists(p):
            continue
        res = D.parse_backbone(p)
        s, keys = D.assign(res)
        if not s:
            continue
        el = D.elements(s, keys)
        for i in keys:
            ss.append({"domain": dom, "pos": i, "ss": el[i][0]})
    return v.merge(pd.DataFrame(ss), on=["domain", "pos"])


def main(rdir, lib3, afdir, out_path):
    ind = pd.read_csv(f"{rdir}/tsuboyama_indels.csv")
    feat = pd.read_csv(f"{rdir}/tsub_features.csv")
    ss = pd.read_csv(f"{rdir}/tsub_ss.csv")
    ss["domain"] = ss["domain"].str.replace(".pdb", "", regex=False)
    wt = json.load(open(f"{rdir}/tsub_wt_dG.json"))
    ind = ind[(ind["kind"] == "del") & (ind["dG_ci"] <= 1.0)].copy()
    ind["wt"] = ind["domain"].map(wt)
    ind = ind[ind["wt"].notna()]
    ind["del"] = ind["wt"] - ind["dG"]
    df = (feat[["domain", "pos", "struct_pos", "plddt", "rsa"]]
          .merge(ss[["domain", "struct_pos", "aa", "ss"]], on=["domain", "struct_pos"])
          .merge(ind.groupby(["domain", "pos"])["del"].mean().reset_index(),
                 on=["domain", "pos"]))
    df["pace"] = df["aa"].map(PACE)
    df["vol"] = df["aa"].map(VOL)
    df = df.dropna(subset=["del", "pace", "vol", "rsa"])

    rep = {"n": int(len(df)), "n_domains": int(df["domain"].nunique()),
           "residues_present": int(df["aa"].nunique()),
           "counts_by_class": {str(k): int(v) for k, v in df["ss"].value_counts().items()}}

    # power check (kill condition 2)
    ok = 0
    for aa, g in df.groupby("aa"):
        if (g["ss"] == "H").sum() >= 100 and (g["ss"] == "L").sum() >= 100:
            ok += 1
    rep["residues_with_100_in_both_H_and_L"] = ok
    rep["D2_underpowered"] = bool(ok < 15)

    # D-1
    import statsmodels.formula.api as smf
    df["dom"] = df["domain"].astype("category")
    cl = {"groups": df["dom"].cat.codes}
    base = smf.ols("Q('del') ~ C(dom) + C(aa) + C(ss) + rsa + vol", data=df).fit(
        cov_type="cluster", cov_kwds=cl)
    full = smf.ols("Q('del') ~ C(dom) + C(aa)*C(ss) + rsa + vol", data=df).fit(
        cov_type="cluster", cov_kwds=cl)
    f = full.compare_f_test(base)
    rep["D1_interaction_F"] = {"F": round(float(f[0]), 3), "p": float(f[1]),
                               "df": float(f[2]), "pass": bool(f[1] < 0.001)}

    # D-2
    hl = contrast_by_residue(df, "H", "L")
    el = contrast_by_residue(df, "E", "L")
    rep["helix_minus_loop_by_residue"] = {k: round(v["mean"], 3) for k, v in hl.items()}
    rep["strand_minus_loop_by_residue"] = {k: round(v["mean"], 3) for k, v in el.items()}
    aas = sorted(hl)
    x = np.array([PACE[a] for a in aas])
    y = np.array([hl[a]["mean"] for a in aas])
    rho, p = stats.spearmanr(x, y)
    rep["D2"] = {"n_residues": len(aas), "spearman_pace_vs_helix_contrast": round(float(rho), 3),
                 "p": float(p), "pass": bool(rho <= -0.4 and p < 0.05)}
    # without Gly and Pro (kill condition 3)
    keep = [a for a in aas if a not in ("G", "P")]
    if len(keep) >= 10:
        r2, p2 = stats.spearmanr([PACE[a] for a in keep], [hl[a]["mean"] for a in keep])
        rep["D2_without_Gly_Pro"] = {"n": len(keep), "spearman": round(float(r2), 3),
                                     "p": float(p2)}
    # volume comparison (D-4)
    rv, pv = stats.spearmanr([VOL[a] for a in aas], y)
    rep["D4_volume_comparison"] = {"spearman_volume_vs_helix_contrast": round(float(rv), 3),
                                   "p": float(pv),
                                   "propensity_beats_volume": bool(abs(rho) > abs(rv))}

    # D-3 within rsa tertiles
    df["t"] = pd.qcut(df["rsa"], 3, labels=["buried", "mid", "exposed"], duplicates="drop")
    strata = {}
    for t, g in df.groupby("t", observed=True):
        h = contrast_by_residue(g, "H", "L")
        a2 = sorted(h)
        if len(a2) >= 8:
            r3, p3 = stats.spearmanr([PACE[a] for a in a2], [h[a]["mean"] for a in a2])
            strata[str(t)] = {"n_residues": len(a2), "spearman": round(float(r3), 3),
                              "p": float(p3)}
    rep["D3_by_burial"] = strata

    # kill condition 1: permute residue identity within domain
    rng = np.random.default_rng(3)
    p = df.copy()
    p["aa"] = p.groupby("domain")["aa"].transform(lambda s: rng.permutation(s.values))
    p["dom"] = p["domain"].astype("category")
    clp = {"groups": p["dom"].cat.codes}
    b2 = smf.ols("Q('del') ~ C(dom) + C(aa) + C(ss) + rsa + vol", data=p).fit(
        cov_type="cluster", cov_kwds=clp)
    f2 = smf.ols("Q('del') ~ C(dom) + C(aa)*C(ss) + rsa + vol", data=p).fit(
        cov_type="cluster", cov_kwds=clp)
    fp = f2.compare_f_test(b2)
    rep["kill1_permutation"] = {"F": round(float(fp[0]), 3), "p": float(fp[1]),
                                "vanishes": bool(fp[1] > 0.01)}

    # D-5 designed replication
    try:
        dd = build_designed(lib3, afdir)
        dd["pace"] = dd["aa"].map(PACE)
        h = contrast_by_residue(dd, "H", "L", nmin=1)
        a3 = sorted(h)
        if len(a3) >= 8:
            r5, p5 = stats.spearmanr([PACE[a] for a in a3], [h[a]["mean"] for a in a3])
            rep["D5_designed"] = {"n_positions": int(len(dd)), "n_residues": len(a3),
                                  "spearman": round(float(r5), 3), "p": float(p5),
                                  "contrasts": {k: round(v["mean"], 3) for k, v in h.items()}}
        else:
            rep["D5_designed"] = {"status": "too few residues with enough domains",
                                  "n_positions": int(len(dd))}
    except Exception as e:  # noqa: BLE001
        rep["D5_designed"] = {"error": str(e)}

    json.dump(rep, open(out_path, "w"), indent=2)
    print(json.dumps(rep, indent=2))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
