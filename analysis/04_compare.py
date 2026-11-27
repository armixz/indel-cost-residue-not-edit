"""Gate 1: ProDomino against the simple baselines on the same measured labels.

Pass marks were fixed in ../RESULTS.md before any score was computed. AraC is excluded from the
primary comparison because the ProDomino authors used it to choose their checkpoint.

Usage: python -I 04_compare.py <baselines_sites.csv> <prodomino.csv> <baselines.json> <out.json>
"""
import json
import sys

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
from scipy.stats import spearmanr

HELD_OUT = ["Flp", "SigF", "TVMV"]


def auroc(y, s):
    m = np.isfinite(s)
    if len(np.unique(y[m])) < 2:
        return None
    return float(roc_auc_score(y[m], s[m]))


def main(sites_path, pd_path, base_path, out_path):
    sites = pd.read_csv(sites_path)
    pdsc = pd.read_csv(pd_path)
    base = json.load(open(base_path))
    rep = {"offset_check": {}, "prodomino_auroc": {}, "comparison": {}}

    # an off-by-one in the residue convention would change everything, so check shifts first
    for off in (-2, -1, 0, 1, 2):
        shifted = pdsc.copy()
        shifted["res0"] = shifted["res0"] - off
        m = sites.merge(shifted[["protein", "res0", "prodomino"]], on=["protein", "res0"], how="left")
        rep["offset_check"][str(off)] = {
            p: (round(a, 3) if (a := auroc(m.loc[m.protein == p, "y"].values,
                                           m.loc[m.protein == p, "prodomino"].values)) else None)
            for p in sorted(sites["protein"].unique())}

    d = sites.merge(pdsc[["protein", "res0", "prodomino"]], on=["protein", "res0"], how="left")
    miss = int(d["prodomino"].isna().sum())
    rep["sites_without_prodomino_score"] = miss

    for p in sorted(d["protein"].unique()):
        g = d[d["protein"] == p]
        a = auroc(g["y"].values, g["prodomino"].values)
        sp = spearmanr(g["prodomino"], g["best_score"], nan_policy="omit")
        rep["prodomino_auroc"][p] = {"auroc": round(a, 3) if a else None,
                                     "spearman_vs_score": round(float(sp.statistic), 3),
                                     "sites": int(len(g)), "positives": int(g["y"].sum())}
    mask = d["protein"].isin(HELD_OUT)
    rep["prodomino_auroc"]["pooled_held_out_3"] = round(auroc(d.loc[mask, "y"].values, d.loc[mask, "prodomino"].values), 3)
    rep["prodomino_auroc"]["pooled_all_4"] = round(auroc(d["y"].values, d["prodomino"].values), 3)

    comb = base["combined_lofo_auroc"]
    singles = base["single_feature_auroc"]
    best_single = {p: max(((f, v[p]) for f, v in singles.items() if v[p]), key=lambda t: t[1])
                   for p in HELD_OUT + ["AraC"]}
    wins_comb = sum(1 for p in HELD_OUT if rep["prodomino_auroc"][p]["auroc"] > comb[p])
    wins_single = sum(1 for p in HELD_OUT if rep["prodomino_auroc"][p]["auroc"] > best_single[p][1])
    rep["comparison"] = {
        "per_protein": {p: {"prodomino": rep["prodomino_auroc"][p]["auroc"],
                            "combined_baseline_lofo": comb[p],
                            "best_single_feature": best_single[p][0],
                            "best_single_auroc": best_single[p][1]}
                        for p in ["AraC"] + HELD_OUT},
        "pooled_held_out_3": {"prodomino": rep["prodomino_auroc"]["pooled_held_out_3"],
                              "combined_baseline_lofo": comb["pooled_held_out_3"]},
        "held_out_proteins_where_prodomino_beats_combined": wins_comb,
        "held_out_proteins_where_prodomino_beats_best_single": wins_single,
    }
    rep["verdict"] = {
        "mark1_beats_combined_on_at_least_2_of_3": bool(wins_comb >= 2),
        "mark2_pooled_held_out_higher": bool(rep["prodomino_auroc"]["pooled_held_out_3"] > comb["pooled_held_out_3"]),
        "passed": bool(wins_comb >= 2 and rep["prodomino_auroc"]["pooled_held_out_3"] > comb["pooled_held_out_3"]),
    }
    rep["arac_reproduction"] = {"published_auroc": 0.84, "measured_here": rep["prodomino_auroc"]["AraC"]["auroc"],
                                "conservation_entropy_auroc": singles["entropy"]["AraC"]}
    d.to_csv(out_path.replace(".json", "_merged.csv"), index=False)
    json.dump(rep, open(out_path, "w"), indent=1)
    print(json.dumps(rep, indent=1))


if __name__ == "__main__":
    main(*sys.argv[1:5])
