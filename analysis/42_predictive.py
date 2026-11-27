"""Does modelling the residue-by-structure interaction improve prediction of indel cost?

If the interaction reported here is real and not already captured, then a predictor that knows
only the structure class, or only the residue, should be beaten by one that knows both together.
Evaluated leave-one-domain-out, so no domain contributes to its own prediction, and scored
per-domain (Spearman within each domain) as well as pooled, because this project has shown the
two are not interchangeable.

Models, all ridge-regularised linear models on the same rows:
  S     structure class only (H/E/L one-hot)
  R     residue identity only (20 one-hot)
  S+R   both, additive
  S*R   both plus their interaction (the claim)
  +geom adds relative solvent accessibility, pLDDT and distance to the chain terminus
  +geom*ss adds the interaction of those with structure class, as a fairness control

Usage: indel_env python 42_predictive.py <lib3_csv> <afdir> <results_dir> <out_json>
"""
import csv
import importlib.util
import json
import os
import re
import sys

import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
sp = importlib.util.spec_from_file_location("d17", os.path.join(HERE, "17_dssp.py"))
D17 = importlib.util.module_from_spec(sp)
sp.loader.exec_module(D17)
PAT = re.compile(r"^(?P<dom>.+?)_(?P<kind>ins|del)(?P<res>[A-Z])(?P<pos>\d+)$")
AAS = list("ACDEFGHIKLMNPQRSTVWY")
SS = ["H", "E", "L"]


def build(lib3, afdir, rdir):
    csv.field_size_limit(10 ** 7)
    wt, rows = {}, []
    with open(lib3, newline="") as fh:
        for r in csv.DictReader(fh):
            nm = r["name"]
            if nm.rsplit("_", 1)[-1].startswith("wt"):
                wt.setdefault(nm.rsplit("_", 1)[0], []).append(float(r["deltaG"]))
                continue
            m = PAT.match(nm)
            if not m:
                continue
            try:
                rows.append({"domain": m["dom"], "kind": m["kind"], "res": m["res"],
                             "pos": int(m["pos"]), "dG": float(r["deltaG"]),
                             "ci": float(r["deltaG_95CI"])})
            except (ValueError, KeyError):
                continue
    wtm = {k: float(np.mean(v)) for k, v in wt.items()}
    v = pd.DataFrame(rows)
    v = v[v["ci"] <= 1.0].copy()
    v["wtdG"] = v["domain"].map(wtm)
    v = v[v["wtdG"].notna()]
    v["cost"] = v["wtdG"] - v["dG"]

    ss = []
    for dom in sorted(v["domain"].unique()):
        p = os.path.join(afdir, dom.replace(".pdb", "") + ".pdb")
        if not os.path.exists(p):
            continue
        res = D17.parse_backbone(p)
        s, keys = D17.assign(res)
        if not s:
            continue
        el = D17.elements(s, keys)
        n = len(keys)
        lo = min(keys)
        for i in keys:
            ss.append({"domain": dom, "pos": i, "aa": res[i]["aa"], "ss": el[i][0],
                       "dist_term": min(i - lo, lo + n - 1 - i)})
    S = pd.DataFrame(ss)
    feat = pd.read_csv(f"{rdir}/tsub_features.csv")[["domain", "pos", "plddt", "rsa"]].copy()
    # the feature table strips the .pdb suffix from domain names; align before merging
    S["_key"] = S["domain"].str.replace(".pdb", "", regex=False)
    feat["_key"] = feat["domain"].astype(str)
    S = S.merge(feat.drop(columns=["domain"]), on=["_key", "pos"], how="left").drop(
        columns=["_key"])
    d = (S.merge(v[v["kind"] == "del"].groupby(["domain", "pos"])["cost"].mean()
                 .rename("y").reset_index(), on=["domain", "pos"]))
    return d.dropna(subset=["y"])


def design(df, spec):
    cols = []
    if "S" in spec:
        cols.append(pd.get_dummies(df["ss"]).reindex(columns=SS, fill_value=0).astype(float))
    if "R" in spec:
        cols.append(pd.get_dummies(df["aa"]).reindex(columns=AAS, fill_value=0).astype(float))
    if "X" in spec:  # interaction
        inter = pd.get_dummies(df["aa"].astype(str) + "_" + df["ss"].astype(str))
        full = [f"{a}_{s}" for a in AAS for s in SS]
        cols.append(inter.reindex(columns=full, fill_value=0).astype(float))
    if "G" in spec:
        g = df[["rsa", "plddt", "dist_term"]].copy()
        g["plddt"] = g["plddt"].fillna(g["plddt"].median())
        g["rsa"] = g["rsa"].fillna(g["rsa"].median())
        g = (g - g.mean()) / (g.std() + 1e-9)
        cols.append(g.astype(float))
    if "GX" in spec:
        g = df[["rsa", "plddt", "dist_term"]].copy()
        g["plddt"] = g["plddt"].fillna(g["plddt"].median())
        g["rsa"] = g["rsa"].fillna(g["rsa"].median())
        g = (g - g.mean()) / (g.std() + 1e-9)
        for s in SS:
            ind = (df["ss"] == s).astype(float).to_numpy()[:, None]
            cols.append(pd.DataFrame(g.to_numpy() * ind, index=df.index,
                                     columns=[f"{c}_{s}" for c in g.columns]).astype(float))
    X = pd.concat(cols, axis=1)
    return X.to_numpy(dtype=float)


def ridge_fit(X, y, lam=1.0):
    X1 = np.hstack([np.ones((len(X), 1)), X])
    A = X1.T @ X1 + lam * np.eye(X1.shape[1])
    A[0, 0] -= lam
    return np.linalg.solve(A, X1.T @ y)


def main(lib3, afdir, rdir, out_path):
    df = build(lib3, afdir, rdir)
    df = df[df["aa"].isin(AAS) & df["ss"].isin(SS)].reset_index(drop=True)
    rep = {"n": int(len(df)), "domains": int(df["domain"].nunique())}

    specs = {"structure_only": {"S"}, "residue_only": {"R"}, "additive_S+R": {"S", "R"},
             "interaction_SxR": {"S", "R", "X"},
             "additive_plus_geometry": {"S", "R", "G"},
             "interaction_plus_geometry": {"S", "R", "X", "G"},
             "additive_plus_geometry_x_structure": {"S", "R", "G", "GX"}}
    doms = df["domain"].unique()
    out = {}
    for name, spec in specs.items():
        X = design(df, spec)
        y = df["y"].to_numpy(dtype=float)
        pred = np.zeros(len(y))
        for dom in doms:
            te = (df["domain"] == dom).to_numpy()
            tr = ~te
            if tr.sum() < 100 or te.sum() == 0:
                pred[te] = y[tr].mean() if tr.sum() else 0.0
                continue
            beta = ridge_fit(X[tr], y[tr])
            pred[te] = np.hstack([np.ones((te.sum(), 1)), X[te]]) @ beta
        per = []
        for dom in doms:
            te = (df["domain"] == dom).to_numpy()
            if te.sum() >= 10 and np.std(y[te]) > 0 and np.std(pred[te]) > 0:
                per.append(float(stats.spearmanr(y[te], pred[te])[0]))
        out[name] = {"n_features": int(X.shape[1]),
                     "rmse": round(float(np.sqrt(np.mean((y - pred) ** 2))), 4),
                     "pooled_spearman": round(float(stats.spearmanr(y, pred)[0]), 4),
                     "per_domain_spearman_mean": round(float(np.mean(per)), 4),
                     "per_domain_spearman_median": round(float(np.median(per)), 4),
                     "per_domain_lower_decile": round(float(np.quantile(per, 0.1)), 4),
                     "n_domains_scored": len(per)}
        out[name]["_per"] = per
    # paired comparison of the two models that matter
    a = np.array(out["additive_S+R"]["_per"])
    b = np.array(out["interaction_SxR"]["_per"])
    n = min(len(a), len(b))
    rng = np.random.default_rng(0)
    diff = b[:n] - a[:n]
    bs = diff[rng.integers(0, n, size=(20000, n))].mean(axis=1)
    rep["paired_gain_interaction_over_additive"] = {
        "mean_delta_per_domain_spearman": round(float(diff.mean()), 4),
        "ci95": [round(float(np.quantile(bs, .025)), 4), round(float(np.quantile(bs, .975)), 4)],
        "domains_improved": int((diff > 0).sum()), "n_domains": n}
    for k in out:
        out[k].pop("_per", None)
    rep["models"] = out
    json.dump(rep, open(out_path, "w"), indent=2)
    print(json.dumps(rep, indent=2))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
