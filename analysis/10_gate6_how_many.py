"""Gate 6: how many measured positions does one protein need before its own model is useful,
and does choosing the panel beat choosing it at random? See ../RESULTS.md.

Also Gate 7 (consensus of ranks) and Gate 8 (ESM-2 features) on the same eight datasets.

Usage: python -I 10_gate6_how_many.py <insertion-lab root> <out.json>
"""
import json
import re
import sys

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.mixture import GaussianMixture
from scipy.stats import rankdata

FEATS = ["plddt", "rsa", "is_loop", "dist_term", "prodomino"]
KS = [8, 16, 24, 32, 48, 64, 96, 128]
REPEATS = 50
PDB3 = dict(zip(["ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", "ILE", "LEU", "LYS",
                 "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"], list("ARNDCEQGHILKMFPSTWYV")))


def kir_map(p):
    """Kir2.1 construct numbering -> AlphaFold residue (FLAG tag DYKDDDDK sits after 115)."""
    return p if p <= 115 else (None if p <= 123 else p - 8)


def auroc(y, s):
    m = np.isfinite(s)
    if m.sum() < 10 or len(np.unique(y[m])) < 2:
        return np.nan
    return float(roc_auc_score(y[m], s[m]))


def build_datasets(root):
    """Returns {name: DataFrame with pos, y, features...}; reuses earlier gate outputs."""
    import importlib.util
    spec = importlib.util.spec_from_file_location("g3", f"{root}/analysis/07_gate3.py")
    g3 = importlib.util.module_from_spec(spec); spec.loader.exec_module(g3)
    md, ds = f"{root}/data/mavedb", f"{root}/data/DI_screen/analysis/input_data"
    pd3 = pd.read_csv(f"{root}/results/prodomino_gate3.csv")
    pdk = pd.read_csv(f"{root}/results/prodomino_kir21.csv")
    pdd = pd.read_csv(f"{root}/results/prodomino.csv")      # DI_screen proteins, res0 0-based
    try:
        esm = pd.read_csv(f"{root}/results/esm_features.csv")
    except Exception:
        esm = None
    out = {}

    def feats_from_pdb(pdb, positions, mapper=None):
        res = g3.parse_pdb(pdb)
        rsa, loop = g3.rsa_per_residue(res), g3.ss_loop(res)
        plddt = {i: float(np.mean(d["b"])) for i, d in res.items()}
        L = max(res)
        rows = []
        for p in positions:
            s = mapper(int(p)) if mapper else int(p)
            rows.append({"pos": int(p), "plddt": plddt.get(s, np.nan) if s else np.nan,
                         "rsa": rsa.get(s, np.nan) if s else np.nan,
                         "is_loop": loop.get(s, np.nan) if s else np.nan,
                         "dist_term": min(int(p) - 1, L - int(p))})
        return pd.DataFrame(rows)

    def cut(x):
        g = GaussianMixture(2, random_state=0, n_init=5).fit(np.asarray(x).reshape(-1, 1))
        return float(np.mean(sorted(g.means_.ravel())))

    # short insertions
    for nm, n in (("Kir2.1_ins1", 1), ("Kir2.1_ins3", 3)):
        t = g3.insertions(f"{md}/kir21_scores.csv")
        t = t[t["pos"].notna()]
        tb = t.assign(nres=1)  # g3.insertions keeps single-residue only; rebuild for n>1
        raw = []
        for r in pd.read_csv(f"{md}/kir21_scores.csv").to_dict("records"):
            m = re.match(r"p\.[A-Za-z]{3}(\d+)_[A-Za-z]{3}\d+ins([A-Za-z]+)$", str(r.get("hgvs_pro") or ""))
            if m and np.isfinite(r.get("score")) and len(re.findall(r"[A-Z][a-z]{2}", m.group(2))) == n:
                raw.append({"pos": int(m.group(1)), "score": float(r["score"])})
        tb = pd.DataFrame(raw).groupby("pos")["score"].median().reset_index()
        f = feats_from_pdb(f"{md}/af_P35561.pdb", tb["pos"], kir_map)
        d = tb.merge(f, on="pos")
        d["prodomino"] = d["pos"].map(pdk.set_index("pos")["prodomino"])
        d["y"] = (d["score"] >= cut(d["score"].values)).astype(int)
        out[nm] = d
    tb = g3.insertions(f"{md}/scores_00000100-a-1.csv").groupby("pos")["score"].median().reset_index()
    f = feats_from_pdb(f"{md}/af_P62593.pdb", tb["pos"])
    d = tb.merge(f, on="pos")
    d["prodomino"] = d["pos"].map(pd3[pd3["name"] == "TEM1"].set_index("pos")["prodomino"])
    d["y"] = (d["score"] >= cut(d["score"].values)).astype(int)
    out["TEM1_ins1"] = d
    dh = pd.concat([g3.insertions(f"{md}/scores_00001219-a-{i}.csv") for i in (1, 2, 3, 4)])
    dh = dh.groupby("pos")["score"].mean().reset_index()
    f = feats_from_pdb(f"{md}/af_P00374.pdb", dh["pos"])
    d = dh.merge(f, on="pos")
    d["prodomino"] = d["pos"].map(pd3[pd3["name"] == "DHFR"].set_index("pos")["prodomino"])
    d["y"] = (d["score"] >= cut(d["score"].values)).astype(int)
    out["DHFR_ins1"] = d
    # domain insertions
    fe = pd.read_csv(f"{root}/results/features.csv")
    for prot, g in fe.groupby("protein"):
        piv = g.pivot_table(index="res0", columns="domain", values="score")
        keep = piv.notna().any(axis=1)
        base = g.drop_duplicates("res0").set_index("res0").reindex(piv.index)
        d = pd.DataFrame({"pos": piv.index[keep] + 1,
                          "y": (piv[keep] > 0).any(axis=1).astype(int).values,
                          "plddt": base.loc[keep, "plddt"].values,
                          "rsa": base.loc[keep, "rsa"].values,
                          "is_loop": (base.loc[keep, "ss"] == "L").astype(float).values,
                          "dist_term": base.loc[keep, "dist_term"].values})
        d["prodomino"] = d["pos"].map(pdd[pdd["protein"] == prot].set_index(pdd[pdd["protein"] == prot]["res0"] + 1)["prodomino"])
        out[f"{prot}_domain"] = d
    if esm is not None:
        key = {"Kir2.1_ins1": "mKir2.1", "Kir2.1_ins3": "mKir2.1", "TEM1_ins1": "TEM1", "DHFR_ins1": "DHFR",
               "AraC_domain": "AraC", "Flp_domain": "Flp", "SigF_domain": "SigF", "TVMV_domain": "TVMV"}
        for nm, d in out.items():
            e = esm[esm["name"] == key[nm]].set_index("pos")
            d["esm_entropy"] = d["pos"].map(e["esm_entropy"])
            d["esm_wt_loglik"] = d["pos"].map(e["esm_wt_loglik"])
    return out


def pick(strategy, X, k, rng):
    n = len(X)
    if strategy == "random":
        return rng.choice(n, min(k, n), replace=False)
    if strategy == "stratified":
        q = pd.qcut(pd.Series(X[:, 0]).rank(method="first"), 10, labels=False)
        idx = []
        per = max(1, k // 10)
        for b in range(10):
            c = np.flatnonzero(q == b)
            if len(c):
                idx += list(rng.choice(c, min(per, len(c)), replace=False))
        idx = list(dict.fromkeys(idx))
        rest = [i for i in range(n) if i not in idx]
        rng.shuffle(rest)
        return np.array((idx + rest)[:k])
    # diverse: greedy max-min distance
    Z = (X - np.nanmean(X, 0)) / (np.nanstd(X, 0) + 1e-9)
    Z = np.nan_to_num(Z)
    sel = [int(rng.integers(n))]
    d = np.linalg.norm(Z - Z[sel[0]], axis=1)
    while len(sel) < min(k, n):
        j = int(np.argmax(d)); sel.append(j)
        d = np.minimum(d, np.linalg.norm(Z - Z[j], axis=1))
    return np.array(sel)


def main(root, out_path):
    data = build_datasets(root)
    rep = {"datasets": {k: {"n": int(len(v)), "positives": int(v["y"].sum())} for k, v in data.items()},
           "zero_measurement": {}, "curves": {}, "gate7": {}, "gate8": {}}

    single = ["plddt", "rsa", "is_loop", "dist_term", "prodomino"]
    if all("esm_entropy" in d for d in data.values()):
        single += ["esm_entropy", "esm_wt_loglik"]
    ORI = {"plddt": -1, "rsa": +1, "is_loop": +1, "dist_term": -1, "prodomino": +1,
           "esm_entropy": +1, "esm_wt_loglik": -1}
    for nm, d in data.items():
        rep["zero_measurement"][nm] = {f: (round(a, 3) if np.isfinite(a := auroc(d["y"].values, ORI[f] * d[f].values.astype(float))) else None)
                                       for f in single if f in d}
        # Gate 7: mean rank of the four robust features
        cols = [c for c in ("plddt", "rsa", "dist_term", "prodomino") if c in d]
        R = np.mean([rankdata(ORI[c] * np.nan_to_num(d[c].values.astype(float), nan=np.nanmedian(d[c]))) for c in cols], 0)
        rep["zero_measurement"][nm]["consensus_rank"] = round(auroc(d["y"].values, R), 3)

    # Gate 6 curves
    for nm, d in data.items():
        X = d[[c for c in FEATS if c in d]].astype(float)
        X = X.fillna(X.median()).values
        y = d["y"].values.astype(int)
        rep["curves"][nm] = {}
        for strat in ("random", "stratified", "diverse"):
            rep["curves"][nm][strat] = {}
            for k in KS:
                if k >= len(d) - 20:
                    continue
                vals = []
                for r in range(REPEATS):
                    rng = np.random.default_rng(1000 * k + r)
                    tr = pick(strat, X, k, rng)
                    te = np.setdiff1d(np.arange(len(d)), tr)
                    if len(np.unique(y[tr])) < 2:
                        continue
                    m = make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000, C=1.0))
                    m.fit(X[tr], y[tr])
                    vals.append(auroc(y[te], m.predict_proba(X[te])[:, 1]))
                if vals:
                    rep["curves"][nm][strat][str(k)] = {"median": round(float(np.nanmedian(vals)), 3),
                                                        "iqr": [round(float(np.nanquantile(vals, .25)), 3),
                                                                round(float(np.nanquantile(vals, .75)), 3)],
                                                        "n_rep": len(vals)}
    # marks
    m1 = []
    for nm in data:
        base = rep["zero_measurement"][nm].get("plddt")
        got = rep["curves"][nm]["random"].get("24", {}).get("median")
        if base is not None and got is not None:
            m1.append((nm, got, base, got > base))
    rep["M1"] = {"per_dataset": [{"dataset": a, "fitted_k24": b, "plddt_only": c, "better": d} for a, b, c, d in m1],
                 "count_better": int(sum(d for *_, d in m1)), "of": len(m1), "passed": bool(sum(d for *_, d in m1) >= 5)}
    diffs = []
    for nm in data:
        r = rep["curves"][nm]["random"].get("24", {}).get("median")
        best = max([rep["curves"][nm][s].get("24", {}).get("median", -1) for s in ("stratified", "diverse")])
        if r is not None and best > 0:
            diffs.append(best - r)
    rep["M2"] = {"mean_gain_at_k24": round(float(np.mean(diffs)), 4) if diffs else None,
                 "per_dataset_gain": [round(x, 3) for x in diffs],
                 "passed": bool(diffs and np.mean(diffs) >= 0.02)}
    json.dump(rep, open(out_path, "w"), indent=1)
    print(json.dumps({k: v for k, v in rep.items() if k != "curves"}, indent=1))
    print("\ncurves (median AUROC, random selection):")
    for nm in rep["curves"]:
        print(" ", nm, {k: v["median"] for k, v in rep["curves"][nm]["random"].items()})


if __name__ == "__main__":
    main(*sys.argv[1:3])
