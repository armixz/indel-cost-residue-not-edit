"""Baseline benchmark: how well do simple features predict measured insertion tolerance?

Orientations are fixed a priori (see ../RESULTS.md): more variable, more exposed, in a loop, lower
pLDDT and nearer a terminus are all expected to be MORE tolerant. AUROC is reported with that
orientation, so a value below 0.5 means the expectation was wrong, and is reported as such.

Usage: python -I 02_baselines.py <features.csv> <out.json>
"""
import json
import sys

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from scipy.stats import spearmanr

FEATURES = {"entropy": +1, "rsa": +1, "is_loop": +1, "plddt": -1, "dist_term": -1}


def build(path):
    df = pd.read_csv(path)
    df["is_loop"] = (df["ss"] == "L").astype(float)
    # AraC: pooled over its five inserted domains, positive if enriched with any of them
    rows = []
    for prot, g in df.groupby("protein"):
        piv = g.pivot_table(index="res0", columns="domain", values="score")
        feat = g.drop_duplicates("res0").set_index("res0").reindex(piv.index)
        y = (piv > 0).any(axis=1) & piv.notna().any(axis=1)
        keep = piv.notna().any(axis=1)
        best = piv.max(axis=1)
        out = feat.loc[keep, list(FEATURES) + ["aa", "length"]].copy()
        out["y"] = y[keep].astype(int)
        out["best_score"] = best[keep]
        out["protein"] = prot
        out["n_domains"] = piv.shape[1]
        rows.append(out.reset_index())
    return pd.concat(rows, ignore_index=True)


def auroc(y, s):
    if len(np.unique(y)) < 2 or np.all(~np.isfinite(s)):
        return None
    m = np.isfinite(s)
    return float(roc_auc_score(y[m], s[m])) if len(np.unique(y[m])) == 2 else None


def main(feat_path, out_path):
    d = build(feat_path)
    prots = sorted(d["protein"].unique())
    rep = {"n_sites": int(len(d)), "by_protein": {}, "single_feature_auroc": {}, "spearman": {}}
    for p in prots:
        g = d[d["protein"] == p]
        rep["by_protein"][p] = {"sites": int(len(g)), "positives": int(g["y"].sum()),
                                "positive_rate": round(float(g["y"].mean()), 3),
                                "domains": int(g["n_domains"].iloc[0])}
    for f, sign in FEATURES.items():
        rep["single_feature_auroc"][f] = {}
        rep["spearman"][f] = {}
        for p in prots:
            g = d[d["protein"] == p]
            a = auroc(g["y"].values, sign * g[f].values)
            rep["single_feature_auroc"][f][p] = round(a, 3) if a else None
            m = np.isfinite(g[f].values) & np.isfinite(g["best_score"].values)
            rep["spearman"][f][p] = round(float(spearmanr(sign * g[f].values[m], g["best_score"].values[m]).statistic), 3) if m.sum() > 5 else None
        vals = [v for v in rep["single_feature_auroc"][f].values() if v]
        rep["single_feature_auroc"][f]["mean"] = round(float(np.mean(vals)), 3) if vals else None

    # combined baseline, leave one protein out
    X = d[list(FEATURES)].astype(float)
    X = X.fillna(X.median())
    y = d["y"].values
    pred = np.full(len(d), np.nan)
    coefs = {}
    for p in prots:
        te = (d["protein"] == p).values
        model = make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000, C=1.0))
        model.fit(X[~te], y[~te])
        pred[te] = model.predict_proba(X[te])[:, 1]
        coefs[p] = {f: round(float(c), 3) for f, c in zip(FEATURES, model[-1].coef_[0])}
    rep["combined_lofo_auroc"] = {p: round(auroc(y[(d["protein"] == p).values], pred[(d["protein"] == p).values]), 3) for p in prots}
    held = [p for p in prots if p != "AraC"]
    mask = d["protein"].isin(held).values
    rep["combined_lofo_auroc"]["pooled_held_out_3"] = round(auroc(y[mask], pred[mask]), 3)
    rep["combined_lofo_auroc"]["pooled_all_4"] = round(auroc(y, pred), 3)
    rep["combined_coefficients_per_fold"] = coefs
    rep["note"] = ("AraC is reported but excluded from the primary comparison because the ProDomino "
                   "authors used it to choose their checkpoint. ProDomino scores are not in this file.")
    d.to_csv(out_path.replace(".json", "_sites.csv"), index=False)
    json.dump(rep, open(out_path, "w"), indent=1)
    print(json.dumps(rep, indent=1))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
