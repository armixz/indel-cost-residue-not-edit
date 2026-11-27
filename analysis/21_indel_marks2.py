"""Marks B, H1, H2, H3 of PREREG_indel.md -- all statistics computed WITHIN domain first.

Replaces 20_indel_marks.py. The pooled estimator mean(middle) - mean(end) over all positions
failed its pre-registered permutation control (-0.169 [-0.272, -0.057] on shuffled data) because
it is confounded by domain composition: domains differ in overall indel cost and in how many
middle versus end positions they contribute, and a within-domain shuffle preserves that. Every
statistic here is therefore computed inside each domain and then averaged over domains, with a
percentile bootstrap over the per-domain vector. Marks are exactly as pre-registered.

cost = dG(wild type) - dG(variant). Larger cost = more destabilizing.

Usage: run.sh 21_indel_marks2.py <results_dir> <out_json>
"""
import json
import sys

import numpy as np
import pandas as pd
from scipy import stats

CI_MAX = 1.0
NBOOT = 10000
RNG = np.random.default_rng(0)


def agg(vals, nboot=NBOOT):
    """Mean of a per-domain vector with a percentile bootstrap CI over domains."""
    v = np.asarray([x for x in vals if x is not None and np.isfinite(x)], dtype=float)
    if len(v) < 10:
        return {"n_domains": int(len(v)), "mean": None, "ci95": [None, None]}
    idx = RNG.integers(0, len(v), size=(nboot, len(v)))
    bs = v[idx].mean(axis=1)
    return {"n_domains": int(len(v)), "mean": round(float(v.mean()), 3),
            "median": round(float(np.median(v)), 3),
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

    df = feat[["domain", "pos", "struct_pos", "plddt", "rsa", "dist_term"]].merge(
        ss[["domain", "struct_pos", "aa", "ss", "idx_end", "idx_n", "idx_c", "elem_len", "d_loop"]],
        on=["domain", "struct_pos"], how="inner")
    for key, kind, res in (("ins_A", "ins", "A"), ("ins_G", "ins", "G"), ("del", "del", None)):
        sub = ind[ind["kind"] == kind]
        if res is not None:
            sub = sub[sub["res"] == res]
        t = sub.groupby(["domain", "pos"])["cost"].mean().rename(key).reset_index()
        df = df.merge(t, on=["domain", "pos"], how="left")
    return df


# ------------------------------------------------------------------ H1
def h1_per_domain(df, col):
    h = df[(df["ss"] == "H") & df[col].notna()]
    out = []
    for _, g in h.groupby("domain"):
        mid, end = g[g["idx_end"] >= 4][col], g[g["idx_end"] <= 1][col]
        if len(mid) >= 1 and len(end) >= 1:
            out.append(float(mid.mean() - end.mean()))
    return out


def mark_H1(df, col, rep, key):
    h = df[(df["ss"] == "H") & df[col].notna()]
    r = agg(h1_per_domain(df, col))
    res = {"n_middle": int((h["idx_end"] >= 4).sum()), "n_end": int((h["idx_end"] <= 1).sum()),
           "mean_cost_middle": round(float(h[h["idx_end"] >= 4][col].mean()), 3),
           "mean_cost_end": round(float(h[h["idx_end"] <= 1][col].mean()), 3),
           "paired_difference_kcal": r}
    d, ci = r["mean"], r["ci95"]
    res["mark_H1a_pass"] = bool(d is not None and d >= 1.0 and ci[0] > 0.5)
    res["mark_H1b_vetter_2.3_in_CI"] = bool(d is not None and ci[0] <= 2.3 <= ci[1])

    strata, allpos = {}, True
    for var in ("rsa", "plddt"):
        s = df[df[var].notna()].copy()
        s["t"] = pd.qcut(s[var], 3, labels=["low", "mid", "high"], duplicates="drop")
        per = {}
        for t, g in s.groupby("t", observed=True):
            a = agg(h1_per_domain(g, col), nboot=4000)
            per[str(t)] = a
            if a["mean"] is None or a["mean"] <= 0:
                allpos = False
        strata[var] = per
    res["H1c_strata"] = strata
    res["mark_H1c_pass"] = bool(allpos)
    res["H1_overall_pass"] = bool(res["mark_H1a_pass"] and allpos)
    rep[key] = res


# ------------------------------------------------------------------ H2
def h2_per_domain(df, col, subset=None):
    h = df[(df["ss"] == "H") & df[col].notna() & df["del"].notna()].copy()
    if subset is not None:
        h = h[subset(h)]
    h["asym"] = h[col] - h["del"]
    rhos, means, near, far, paired = [], [], [], [], []
    for _, g in h.groupby("domain"):
        if len(g) >= 10 and g["d_loop"].std() > 0 and g["asym"].std() > 0:
            rhos.append(float(stats.spearmanr(g["d_loop"], g["asym"])[0]))
        means.append(float(g["asym"].mean()))
        n, f = g[g["d_loop"] <= 2]["asym"], g[g["d_loop"] >= 5]["asym"]
        if len(n):
            near.append(float(n.mean()))
        if len(f):
            far.append(float(f.mean()))
        if len(n) and len(f):
            paired.append(float(f.mean() - n.mean()))
    return rhos, means, near, far, paired


def mark_H2(df, col, rep, key):
    rhos, means, near, far, paired = h2_per_domain(df, col)
    r = agg(rhos)
    res = {"spearman_asym_vs_dloop_per_domain": r,
           "mean_asymmetry_ins_minus_del": agg(means),
           "near_loop_le2": agg(near), "far_from_loop_ge5": agg(far),
           "paired_far_minus_near": agg(paired)}
    d, ci = r["mean"], r["ci95"]
    res["mark_H2a_pass"] = bool(d is not None and abs(d) >= 0.15 and (ci[0] > 0 or ci[1] < 0))
    p = res["paired_far_minus_near"]
    res["mark_H2b_pass"] = bool(res["near_loop_le2"]["mean"] is not None
                                and res["near_loop_le2"]["mean"] < 0
                                and res["far_from_loop_ge5"]["mean"] > 0
                                and p["ci95"][0] is not None and p["ci95"][0] > 0)

    ar, am, an, af, ap = h2_per_domain(df, col, subset=lambda h: h["aa"] == "A")
    nA = sum(len(x) for x in [am])
    if len(am) < 30 or nA < 30:
        res["H2c"] = {"status": "UNTESTED - alanine subset too small", "n_domains": len(am)}
        res["mark_H2c_pass"] = None
        res["H2_overall"] = "untested"
    else:
        aR = agg(ar)
        res["H2c_alanine_only"] = {"spearman": aR, "mean_asymmetry": agg(am)}
        same = (aR["mean"] is not None and d is not None and np.sign(aR["mean"]) == np.sign(d))
        res["mark_H2c_pass"] = bool(same and (aR["ci95"][0] > 0 or aR["ci95"][1] < 0))
        res["H2_overall"] = "pass" if (res["mark_H2a_pass"] and res["mark_H2c_pass"]) else "fail"
    rep[key] = res


# ------------------------------------------------------------------ B
def mark_B(df, rep):
    out = {}
    for target in ("ins_A", "ins_G", "del"):
        sub = df[df[target].notna()]
        res = {}
        for pred in ("plddt", "rsa", "dist_term"):
            s = sub[sub[pred].notna()]
            pooled = abs(float(stats.pearsonr(s[pred], s[target])[0]))
            per = [abs(float(stats.pearsonr(g[pred], g[target])[0]))
                   for _, g in s.groupby("domain")
                   if len(g) >= 10 and g[pred].std() > 0 and g[target].std() > 0]
            per = np.array(per)
            res[pred] = {"pooled_abs_r": round(pooled, 3),
                         "per_domain_median_abs_r": round(float(np.median(per)), 3),
                         "per_domain_lower_decile": round(float(np.quantile(per, 0.1)), 3),
                         "n_domains": int(len(per)),
                         "pooled_minus_per_domain": round(pooled - float(np.median(per)), 3)}
        out[target] = res
    gaps = [v["pooled_minus_per_domain"] for t in out.values() for v in t.values()]
    out["median_gap"] = round(float(np.median(gaps)), 3)
    out["mark_B1_pass_pooled_higher_by_0.10"] = bool(np.median(gaps) >= 0.10)
    out["note"] = ("mark B1 as written asked whether pooling INFLATES r. For these single "
                   "features the gap is negative: pooling across domains adds between-domain "
                   "variance in cost that a within-domain feature cannot explain, so pooled r is "
                   "LOWER. A model that also predicts domain-level stability would gain from the "
                   "same variance. The direction of the bias therefore depends on the predictor, "
                   "which makes comparing a pooled number with a per-domain one worse than simply "
                   "optimistic.")
    rep["B_metric_audit"] = out


# ------------------------------------------------------------------ H3
def mark_H3(df, rep):
    h = df[df["ins_A"].notna() & df["del"].notna()].copy()
    h["asym"] = h["ins_A"] - h["del"]
    out = {}
    for s, g in h.groupby("ss"):
        out[str(s)] = agg([float(x["asym"].mean()) for _, x in g.groupby("domain")])
    rep["H3_exploratory"] = {"note": "declared exploratory in advance; positive = insertion worse",
                             "by_ss": out}


def permutation(df, rep, seed=1):
    rng = np.random.default_rng(seed)
    p = df.copy()
    for col in ("ins_A", "del"):
        p[col] = p.groupby("domain")[col].transform(lambda s: rng.permutation(s.values))
    sub = {}
    mark_H1(p, "ins_A", sub, "H1p")
    mark_H2(p, "ins_A", sub, "H2p")
    h1, h2 = sub["H1p"]["paired_difference_kcal"], sub["H2p"]["spearman_asym_vs_dloop_per_domain"]
    rep["permutation_control"] = {
        "H1_shuffled": h1, "H1_vanishes": bool(h1["ci95"][0] <= 0 <= h1["ci95"][1]),
        "H2_shuffled": h2, "H2_vanishes": bool(h2["ci95"][0] <= 0 <= h2["ci95"][1])}


def main(rdir, out_path):
    df = load(rdir)
    rep = {"n_positions": int(len(df)), "n_domains": int(df["domain"].nunique()),
           "n_ins_A": int(df["ins_A"].notna().sum()), "n_ins_G": int(df["ins_G"].notna().sum()),
           "n_del": int(df["del"].notna().sum())}
    mark_B(df, rep)
    mark_H1(df, "ins_A", rep, "H1_ala")
    mark_H2(df, "ins_A", rep, "H2_ala")
    permutation(df, rep)
    mark_H3(df, rep)
    mark_H1(df, "ins_G", rep, "H1_gly_replication")
    mark_H2(df, "ins_G", rep, "H2_gly_replication")
    json.dump(rep, open(out_path, "w"), indent=2)
    print(json.dumps(rep, indent=2))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
