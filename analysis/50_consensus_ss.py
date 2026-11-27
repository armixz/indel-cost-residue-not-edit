"""Does the result survive if the loop class is cleaned against experimental structures?

Script 49 compared our secondary-structure assignment against the HELIX/SHEET records of
the deposited experimental entries and found a systematic, one-sided disagreement: in 3,866
positions we say loop and the depositor says helix or strand, against only 619 the other
way. Our helix and strand calls are clean (96% and 90% agree with experiment) but our loop
class contains a few thousand residues that the experimental structure calls structured.

That contamination can only dilute a helix-minus-loop contrast, so the reported effect
should be a lower bound. "Should be" is a prediction, not a result, and this script tests
it instead of asserting it.

Two restricted analyses, both on the 410 domains with usable deposited records:

  consensus  keep only positions where our assignment and the deposited records agree.
             This removes the contamination from both directions at once.
  deposited  use the deposited assignment itself wherever it exists, ignoring ours.
             A harsher test: it replaces our labels rather than filtering them.

If the contrast grows under both, the dilution argument is confirmed and the manuscript's
number is conservative. If it shrinks, the manuscript is overstating and must be corrected.
The direction is declared here before the numbers are read.

Usage: python 50_consensus_ss.py <lib3_csv> <afdir> <results_dir> <out_json>
"""
import importlib.util
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
sp = importlib.util.spec_from_file_location("d43", os.path.join(HERE, "43_discriminate.py"))
D43 = importlib.util.module_from_spec(sp)
sp.loader.exec_module(D43)

RNG = np.random.default_rng(0)


def offset_and_rho(d):
    """Helix-minus-loop deletion-cost contrast per residue, paired within domain, plus
    the correlation of those contrasts with helix propensity. Same estimator as the
    manuscript."""
    cell = {}
    for (dom, aa), g in d.groupby(["domain", "aa"]):
        h, lp = g[g["ss"] == "H"]["cost"], g[g["ss"] == "L"]["cost"]
        if len(h) >= 1 and len(lp) >= 1:
            cell.setdefault(str(aa), []).append(float(h.mean() - lp.mean()))
    pts = {a: (float(np.mean(v)), len(v)) for a, v in cell.items() if len(v) >= 5}
    aas = sorted(a for a in pts if a in D43.PACE)
    if len(aas) < 8:
        return None
    y = np.array([pts[a][0] for a in aas])
    x = np.array([D43.PACE[a] for a in aas])
    rho = float(stats.spearmanr(x, y)[0])
    boot = []
    doms = sorted(d["domain"].unique())
    for _ in range(1000):
        pick = set(RNG.choice(doms, len(doms), replace=True))
        sub = d[d["domain"].isin(pick)]
        c2 = {}
        for (dom, aa), g in sub.groupby(["domain", "aa"]):
            h, lp = g[g["ss"] == "H"]["cost"], g[g["ss"] == "L"]["cost"]
            if len(h) >= 1 and len(lp) >= 1:
                c2.setdefault(str(aa), []).append(float(h.mean() - lp.mean()))
        p2 = {a: float(np.mean(v)) for a, v in c2.items() if len(v) >= 5 and a in D43.PACE}
        if len(p2) >= 8:
            k = sorted(p2)
            r = stats.spearmanr([D43.PACE[a] for a in k], [p2[a] for a in k])[0]
            if np.isfinite(r):
                boot.append(float(r))
    return {"n_residues": len(aas), "spearman": round(rho, 4),
            "spearman_ci95": [round(float(q), 4) for q in np.quantile(boot, [.025, .975])]
            if boot else [None, None],
            "mean_contrast_kcal": round(float(np.mean(y)), 4),
            "n_rows": int(len(d)), "n_domains": int(d["domain"].nunique())}


def main(lib3, afdir, rdir, out_path):
    d = D43.load(lib3, afdir, rdir)
    X = pd.read_csv(os.path.join(rdir, "pdb_crosscheck_residues.csv"))

    rep = {"prediction_declared_before_reading_numbers":
           "cleaning the loop class should INCREASE the magnitude of the helix-minus-loop "
           "contrast, because our loop class is contaminated with residues the experimental "
           "structure calls helix or strand, which dilutes the difference"}

    m = d.merge(X, on=["domain", "pos"], how="inner")
    rep["coverage"] = {
        "rows_in_full_analysis": int(len(d)),
        "rows_with_deposited_records": int(len(m)),
        "domains_with_deposited_records": int(m["domain"].nunique())}

    # (a) the manuscript's own estimator, restricted to the domains that have records,
    #     so the comparison is like for like and not a change of sample
    rep["A_ours_on_covered_domains"] = offset_and_rho(m)

    # (b) consensus positions only
    cons = m[m["ours"] == m["deposited"]]
    rep["B_consensus_positions"] = offset_and_rho(cons)

    # (c) deposited labels substituted for ours
    dep = m.copy()
    dep["ss"] = dep["deposited"]
    rep["C_deposited_labels"] = offset_and_rho(dep)

    base = rep["A_ours_on_covered_domains"]
    out = []
    for k in ("B_consensus_positions", "C_deposited_labels"):
        v = rep[k]
        if base and v:
            out.append({
                "variant": k,
                "contrast_change": round(v["mean_contrast_kcal"]
                                         - base["mean_contrast_kcal"], 4),
                "spearman_change": round(v["spearman"] - base["spearman"], 4)})
    rep["comparison_to_baseline"] = out
    if base and rep["B_consensus_positions"] and rep["C_deposited_labels"]:
        grew = all(abs(rep[k]["mean_contrast_kcal"]) >= abs(base["mean_contrast_kcal"])
                   for k in ("B_consensus_positions", "C_deposited_labels"))
        rep["verdict"] = ("CONFIRMED: cleaning the loop class does not shrink the contrast, "
                          "so the reported effect is conservative" if grew else
                          "NOT CONFIRMED: the contrast shrinks when the loop class is "
                          "cleaned; the manuscript number is not a lower bound")

    json.dump(rep, open(out_path, "w"), indent=2)
    print(json.dumps(rep, indent=2))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
