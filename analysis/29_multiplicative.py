"""Is the insertion/deletion asymmetry just a scale effect?

Across classes the ratio deletion-cost / insertion-cost is 1.29 (strand), 1.15 (helix), 1.14
(loop). If deletions were simply a fixed multiple of insertions, then any position or class with a
higher absolute cost would show a larger ABSOLUTE asymmetry with no structural mechanism at all,
and every result in this project so far would be a restatement of "strand positions are costly".

This script separates the two:
  1. the deletion-versus-insertion relationship fitted WITHIN domain, with and without class terms
  2. asymmetry against insertion cost, binned, per class: does the strand line sit below the helix
     line at the SAME insertion cost?
  3. a log-ratio formulation, which is the right scale if the effect is multiplicative
  4. the same comparison restricted to the 84 mixed alpha/beta domains, where the class contrast
     is within-protein by construction

Usage: indel_env python 29_multiplicative.py <results_dir> <out_json>
"""
import json
import sys

import numpy as np
import pandas as pd
from scipy import stats

CI_MAX = 1.0
RNG = np.random.default_rng(0)


def agg(vals, nboot=20000):
    v = np.asarray([x for x in vals if x is not None and np.isfinite(x)], dtype=float)
    if len(v) < 5:
        return {"n_domains": int(len(v)), "mean": None, "ci95": [None, None]}
    idx = RNG.integers(0, len(v), size=(nboot, len(v)))
    bs = v[idx].mean(axis=1)
    return {"n_domains": int(len(v)), "mean": round(float(v.mean()), 3),
            "ci95": [round(float(np.quantile(bs, 0.025)), 3),
                     round(float(np.quantile(bs, 0.975)), 3)]}


def load(rdir):
    ind = pd.read_csv(f"{rdir}/tsuboyama_indels.csv")
    feat = pd.read_csv(f"{rdir}/tsub_features.csv")
    ss = pd.read_csv(f"{rdir}/tsub_ss.csv")
    ss["domain"] = ss["domain"].str.replace(".pdb", "", regex=False)
    wt = json.load(open(f"{rdir}/tsub_wt_dG.json"))
    ind = ind[ind["dG_ci"] <= CI_MAX].copy()
    ind["wt"] = ind["domain"].map(wt)
    ind = ind[ind["wt"].notna()]
    ind["cost"] = ind["wt"] - ind["dG"]
    df = feat[["domain", "pos", "struct_pos", "plddt", "rsa"]].merge(
        ss[["domain", "struct_pos", "aa", "ss", "idx_end", "elem_len"]],
        on=["domain", "struct_pos"], how="inner")
    for key, kind, res in (("ins_A", "ins", "A"), ("del", "del", None)):
        sub = ind[ind["kind"] == kind]
        if res is not None:
            sub = sub[sub["res"] == res]
        t = sub.groupby(["domain", "pos"])["cost"].mean().rename(key).reset_index()
        df = df.merge(t, on=["domain", "pos"], how="left")
    df = df[df["ins_A"].notna() & df["del"].notna()].copy()
    df["asym"] = df["ins_A"] - df["del"]
    return df


def main(rdir, out_path):
    df = load(rdir)
    rep = {"n_positions": int(len(df))}

    comp = df.groupby("domain")["ss"].value_counts(normalize=True).unstack(fill_value=0.0)
    for c in ("H", "E", "L"):
        if c not in comp:
            comp[c] = 0.0
    mixed = set(comp[(comp["H"] > 0) & (comp["E"] > 0)].index)

    # 1. slope of deletion on insertion, within domain
    slopes = []
    for _, g in df.groupby("domain"):
        if len(g) >= 15 and g["ins_A"].std() > 0:
            slopes.append(float(np.polyfit(g["ins_A"], g["del"], 1)[0]))
    rep["within_domain_slope_del_on_ins"] = agg(slopes)
    rep["pooled_ratio_by_class"] = {
        str(s): {"ins": round(float(g["ins_A"].mean()), 3),
                 "del": round(float(g["del"].mean()), 3),
                 "ratio": round(float(g["del"].mean() / g["ins_A"].mean()), 3)}
        for s, g in df.groupby("ss")}

    # 2. asymmetry against insertion cost, matched bins
    d = df[(df["ins_A"] > 0.05)].copy()
    d["bin"] = pd.qcut(d["ins_A"], 6, duplicates="drop")
    binned = {}
    for b, g in d.groupby("bin", observed=True):
        row = {"ins_mid": round(float(g["ins_A"].mean()), 2), "n": int(len(g))}
        for s in ("E", "H", "L"):
            x = g[g["ss"] == s]
            row[s] = {"n": int(len(x)),
                      "asym": round(float(x["asym"].mean()), 3) if len(x) >= 20 else None}
        binned[str(b)] = row
    rep["asym_by_insertion_cost_bin"] = binned
    # strand minus helix within the same insertion-cost bin, averaged over bins
    diffs = [b["E"]["asym"] - b["H"]["asym"] for b in binned.values()
             if b["E"]["asym"] is not None and b["H"]["asym"] is not None]
    rep["strand_minus_helix_at_matched_insertion_cost"] = {
        "n_bins": len(diffs),
        "mean": round(float(np.mean(diffs)), 3) if diffs else None,
        "per_bin": [round(x, 3) for x in diffs]}

    # 3. log ratio scale
    l = df[(df["ins_A"] > 0.2) & (df["del"] > 0.2)].copy()
    l["logratio"] = np.log(l["del"] / l["ins_A"])
    rep["log_ratio_by_class"] = {
        str(s): agg([float(x["logratio"].mean()) for _, x in g.groupby("domain")
                     if len(x) >= 3])
        for s, g in l.groupby("ss")}
    pl = []
    for _, g in l.groupby("domain"):
        e, h = g[g["ss"] == "E"]["logratio"], g[g["ss"] == "H"]["logratio"]
        if len(e) >= 2 and len(h) >= 2:
            pl.append(float(e.mean() - h.mean()))
    rep["log_ratio_E_minus_H_paired"] = agg(pl)
    pl2 = []
    for _, g in l.groupby("domain"):
        e, o = g[g["ss"] == "E"]["logratio"], g[g["ss"] == "L"]["logratio"]
        if len(e) >= 2 and len(o) >= 2:
            pl2.append(float(e.mean() - o.mean()))
    rep["log_ratio_E_minus_L_paired"] = agg(pl2)

    # 4. mixed domains only, insertion cost relative to the domain's loops
    mm = df[df["domain"].isin(mixed)]
    rel = {"E": [], "H": []}
    for _, g in mm.groupby("domain"):
        base = g[g["ss"] == "L"]["ins_A"]
        if len(base) < 2:
            continue
        b = float(base.mean())
        for c in ("E", "H"):
            x = g[g["ss"] == c]["ins_A"]
            if len(x) >= 2:
                rel[c].append(float(x.mean()) - b)
    rep["mixed_only_insertion_cost_rel_loop"] = {c: agg(v) for c, v in rel.items()}
    pe = []
    for _, g in mm.groupby("domain"):
        e, h = g[g["ss"] == "E"]["ins_A"], g[g["ss"] == "H"]["ins_A"]
        if len(e) >= 2 and len(h) >= 2:
            pe.append(float(e.mean() - h.mean()))
    rep["mixed_only_insertion_cost_E_minus_H_paired"] = agg(pe)
    pd_ = []
    for _, g in mm.groupby("domain"):
        e, h = g[g["ss"] == "E"]["del"], g[g["ss"] == "H"]["del"]
        if len(e) >= 2 and len(h) >= 2:
            pd_.append(float(e.mean() - h.mean()))
    rep["mixed_only_deletion_cost_E_minus_H_paired"] = agg(pd_)

    json.dump(rep, open(out_path, "w"), indent=2)
    print(json.dumps(rep, indent=2))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
