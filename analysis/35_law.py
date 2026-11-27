"""Is there a simple law relating deletion cost to insertion cost at the same position?

29_multiplicative.py found a within-domain slope of deletion on insertion of 0.695
[0.667, 0.723] -- below one -- and that the asymmetry is largest where the insertion is cheapest
(-0.758 in the cheapest sixth of positions, -0.265 in the dearest). Together those imply an
additive offset rather than a multiplicative one: a deletion may carry a roughly fixed extra
penalty that matters most where the position is otherwise tolerant.

Fits, all within domain and summarised over domains:
  1. del = a + b * ins, reporting a and b separately
  2. the same on the designed miniproteins, as an independent set
  3. how well a one-parameter model (del = ins + a) does against the two-parameter fit
  4. the offset by secondary structure class -- this is where the strand result should appear if
     the strand effect is an offset effect rather than a slope effect
  5. a practical statement: at positions in the most tolerant quartile, how much better is an
     insertion than a deletion?

Usage: indel_env python 35_law.py <lib3_csv> <results_dir> <afdir> <out_json>
"""
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
exec(open(os.path.join(HERE, "34_scale.py")).read().split("def main")[0])  # noqa: S102


def fit_per_domain(df, xcol, ycol, nmin=12):
    a_s, b_s, r_s = [], [], []
    for _, g in df.groupby("domain"):
        g = g[g[xcol].notna() & g[ycol].notna()]
        if len(g) < nmin or g[xcol].std() <= 0:
            continue
        b, a = np.polyfit(g[xcol], g[ycol], 1)
        a_s.append(float(a))
        b_s.append(float(b))
        r_s.append(float(np.corrcoef(g[xcol], g[ycol])[0, 1]))
    return a_s, b_s, r_s


def main(lib3, rdir, afdir, out_path):
    D, wtD = designed(lib3, afdir)
    N, wtN = natural(rdir)
    D = D.rename(columns={"ins": "ins"})
    rep = {}

    for name, df in (("natural", N), ("designed", D)):
        a, b, r = fit_per_domain(df, "ins", "del")
        rep[name] = {"intercept_kcal": agg(a), "slope": agg(b), "within_domain_r": agg(r)}
        # one-parameter model: del - ins is a constant offset
        off = [float((g["del"] - g["ins"]).mean()) for _, g in df.groupby("domain")
               if len(g) >= 12]
        rep[name]["mean_offset_del_minus_ins"] = agg(off)
        # residual spread of each model
        rss2, rss1 = [], []
        for _, g in df.groupby("domain"):
            g = g[g["ins"].notna() & g["del"].notna()]
            if len(g) < 12 or g["ins"].std() <= 0:
                continue
            bb, aa_ = np.polyfit(g["ins"], g["del"], 1)
            rss2.append(float(np.mean((g["del"] - (aa_ + bb * g["ins"])) ** 2)))
            c = float((g["del"] - g["ins"]).mean())
            rss1.append(float(np.mean((g["del"] - (g["ins"] + c)) ** 2)))
        rep[name]["rmse_two_param"] = round(float(np.sqrt(np.mean(rss2))), 3)
        rep[name]["rmse_one_param_offset_only"] = round(float(np.sqrt(np.mean(rss1))), 3)

        # 4. offset by class
        rep[name]["offset_by_class"] = {
            str(s): agg([float((x["del"] - x["ins"]).mean()) for _, x in g.groupby("domain")
                         if len(x) >= 3])
            for s, g in df.groupby("ss")}

        # 5. practical: tolerant quartile
        q = df["ins"].quantile(0.25)
        tol = df[df["ins"] <= q]
        rep[name]["tolerant_quartile"] = {
            "insertion_cost_cutoff": round(float(q), 3),
            "n": int(len(tol)),
            "mean_insertion_cost": round(float(tol["ins"].mean()), 3),
            "mean_deletion_cost": round(float(tol["del"].mean()), 3),
            "offset": agg([float((x["del"] - x["ins"]).mean()) for _, x in tol.groupby("domain")
                           if len(x) >= 3])}
        q3 = df["ins"].quantile(0.75)
        intol = df[df["ins"] >= q3]
        rep[name]["intolerant_quartile"] = {
            "mean_insertion_cost": round(float(intol["ins"].mean()), 3),
            "mean_deletion_cost": round(float(intol["del"].mean()), 3),
            "offset": agg([float((x["del"] - x["ins"]).mean()) for _, x in intol.groupby("domain")
                           if len(x) >= 3])}

    json.dump(rep, open(out_path, "w"), indent=2)
    print(json.dumps(rep, indent=2))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
