"""Marks B, H1, H2, H3 of PREREG_indel.md, in the pre-registered order.

cost = -(ddG) = dG(wild type) - dG(variant). Larger cost = more destabilizing.
Every confidence interval is bootstrapped over DOMAINS, never over positions.

Usage: run.sh 20_indel_marks.py <results_dir> <out_json>
"""
import json
import sys

import numpy as np
import pandas as pd
from scipy import stats

CI_MAX = 1.0
NBOOT = 2000


def boot_by_domain(df, fn, n=NBOOT, seed=0):
    rng = np.random.default_rng(seed)
    doms = df["domain"].unique()
    idx = {d: df.index[df["domain"] == d].to_numpy() for d in doms}
    vals = []
    for _ in range(n):
        pick = rng.choice(doms, len(doms), replace=True)
        sub = df.loc[np.concatenate([idx[d] for d in pick])]
        v = fn(sub)
        if v is not None and np.isfinite(v):
            vals.append(v)
    if not vals:
        return [float("nan"), float("nan")]
    return [round(float(q), 3) for q in np.nanquantile(vals, [0.025, 0.975])]


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

    tab = {}
    for kind, res in (("ins", "A"), ("ins", "G"), ("del", None)):
        sub = ind[ind["kind"] == kind]
        if res is not None:
            sub = sub[sub["res"] == res]
        key = f"{kind}_{res}" if res else kind
        tab[key] = sub.groupby(["domain", "pos"])["cost"].mean().rename(key).reset_index()

    base = feat[["domain", "pos", "struct_pos", "plddt", "rsa", "dist_term"]]
    df = base.merge(ss[["domain", "struct_pos", "aa", "ss", "idx_end", "elem_len", "d_loop"]],
                    on=["domain", "struct_pos"], how="inner")
    for k, t in tab.items():
        df = df.merge(t, on=["domain", "pos"], how="left")
    return df


# ---------------------------------------------------------------- B
def mark_B(df, rep):
    out = {}
    for target in ("ins_A", "del"):
        sub = df[df[target].notna()]
        res = {}
        for pred, sign in (("plddt", -1), ("rsa", 1), ("dist_term", -1)):
            s = sub[sub[pred].notna()]
            pooled = float(stats.pearsonr(s[pred] * sign, s[target])[0])
            per = []
            for d, g in s.groupby("domain"):
                if len(g) >= 10 and g[pred].std() > 0 and g[target].std() > 0:
                    per.append(float(stats.pearsonr(g[pred] * sign, g[target])[0]))
            per = np.array(per)
            res[pred] = {"pooled_r": round(pooled, 3),
                         "per_domain_median_r": round(float(np.median(per)), 3),
                         "per_domain_lower_decile": round(float(np.quantile(per, 0.1)), 3),
                         "per_domain_iqr": [round(float(np.quantile(per, 0.25)), 3),
                                            round(float(np.quantile(per, 0.75)), 3)],
                         "n_domains": int(len(per)),
                         "gap_pooled_minus_median": round(pooled - float(np.median(per)), 3)}
        out[target] = res
    gaps = [v["gap_pooled_minus_median"] for t in out.values() for v in t.values()]
    out["mark_B1_pass"] = bool(np.median(gaps) >= 0.10)
    out["median_gap_all"] = round(float(np.median(gaps)), 3)
    rep["B_metric_audit"] = out


# ---------------------------------------------------------------- H1
def h1_diff(s, col="ins_A"):
    """Paired WITHIN-domain estimator.

    The pooled estimator mean(middle) - mean(end) over all positions is confounded by domain
    composition: domains differ in overall indel cost and in how many middle vs end positions
    they contribute, and a within-domain permutation preserves that, so the pre-registered
    permutation control did not vanish for it (it gave -0.169 [-0.272, -0.057]). Computing the
    difference inside each domain first and then averaging over domains removes the composition
    term. The marks of PREREG_indel.md are unchanged; only the estimator is corrected, as the
    pre-registered control requires.
    """
    h = s[(s["ss"] == "H") & s[col].notna()]
    per = []
    for _, g in h.groupby("domain"):
        mid, end = g[g["idx_end"] >= 4][col], g[g["idx_end"] <= 1][col]
        if len(mid) >= 1 and len(end) >= 1:
            per.append(float(mid.mean() - end.mean()))
    if len(per) < 10:
        return np.nan
    return float(np.mean(per))


def mark_H1(df, rep, col="ins_A", key="H1"):
    h = df[(df["ss"] == "H") & df[col].notna()].copy()
    diff = h1_diff(h, col)
    ci = boot_by_domain(h, lambda s: h1_diff(s, col))
    res = {"n_middle": int((h["idx_end"] >= 4).sum()), "n_end": int((h["idx_end"] <= 1).sum()),
           "n_domains": int(h["domain"].nunique()),
           "mean_cost_middle": round(float(h[h["idx_end"] >= 4][col].mean()), 3),
           "mean_cost_end": round(float(h[h["idx_end"] <= 1][col].mean()), 3),
           "difference_kcal": round(diff, 3), "ci95_by_domain": ci,
           "mark_H1a_pass": bool(diff >= 1.0 and ci[0] > 0.5),
           "mark_H1b_vetter_2.3_in_CI": bool(ci[0] <= 2.3 <= ci[1])}

    strata = {}
    allpos = True
    for name, var in (("rsa", "rsa"), ("plddt", "plddt")):
        hs = h[h[var].notna()].copy()
        hs["t"] = pd.qcut(hs[var], 3, labels=["low", "mid", "high"], duplicates="drop")
        per = {}
        for t, g in hs.groupby("t", observed=True):
            d = h1_diff(g, col)
            per[str(t)] = {"n": int(len(g)), "difference_kcal": None if np.isnan(d) else round(d, 3),
                           "ci95_by_domain": boot_by_domain(g.reset_index(drop=True),
                                                            lambda s: h1_diff(s, col), n=600)}
            if np.isnan(d) or d <= 0:
                allpos = False
        strata[name] = per
    res["H1c_strata"] = strata
    res["mark_H1c_pass"] = bool(allpos)
    res["H1_overall_pass"] = bool(res["mark_H1a_pass"] and allpos)
    rep[key] = res


# ---------------------------------------------------------------- H2
def h2_rho(s):
    g = s[s["asym"].notna() & s["d_loop"].notna()]
    if len(g) < 30 or g["d_loop"].std() == 0:
        return np.nan
    return float(stats.spearmanr(g["d_loop"], g["asym"])[0])


def mark_H2(df, rep, col="ins_A", key="H2"):
    h = df[(df["ss"] == "H") & df[col].notna() & df["del"].notna()].copy()
    h["asym"] = h[col] - h["del"]
    rho = h2_rho(h)
    ci = boot_by_domain(h.reset_index(drop=True), h2_rho)
    dmean = lambda s: float(s.groupby("domain")["asym"].mean().mean())  # within-domain, then over domains
    near, far = h[h["d_loop"] <= 2]["asym"], h[h["d_loop"] >= 5]["asym"]
    cin = boot_by_domain(h[h["d_loop"] <= 2].reset_index(drop=True), dmean)
    cif = boot_by_domain(h[h["d_loop"] >= 5].reset_index(drop=True), dmean)
    res = {"n_positions": int(len(h)), "n_domains": int(h["domain"].nunique()),
           "spearman_rho_asym_vs_dloop": None if np.isnan(rho) else round(rho, 3),
           "ci95_by_domain": ci,
           "mark_H2a_pass": bool((not np.isnan(rho)) and abs(rho) >= 0.15
                                 and (ci[0] > 0 or ci[1] < 0)),
           "near_loop_le2": {"n": int(len(near)), "mean_asym": round(float(near.mean()), 3),
                             "ci95": cin},
           "far_from_loop_ge5": {"n": int(len(far)), "mean_asym": round(float(far.mean()), 3),
                                 "ci95": cif}}
    res["mark_H2b_pass"] = bool(near.mean() < 0 < far.mean() and cin[1] < cif[0])

    ala = h[h["aa"] == "A"]
    if ala["domain"].nunique() < 30 or len(ala) < 200:
        res["H2c"] = {"n": int(len(ala)), "n_domains": int(ala["domain"].nunique()),
                      "status": "UNTESTED - subset too small, H2 reported as untested"}
        res["mark_H2c_pass"] = None
        res["H2_overall"] = "untested"
    else:
        ra = h2_rho(ala)
        cia = boot_by_domain(ala.reset_index(drop=True), h2_rho)
        res["H2c"] = {"n": int(len(ala)), "n_domains": int(ala["domain"].nunique()),
                      "spearman_rho": None if np.isnan(ra) else round(ra, 3),
                      "ci95_by_domain": cia}
        same = (not np.isnan(ra)) and (not np.isnan(rho)) and np.sign(ra) == np.sign(rho)
        res["mark_H2c_pass"] = bool(same and (cia[0] > 0 or cia[1] < 0))
        res["H2_overall"] = "pass" if (res["mark_H2a_pass"] and res["mark_H2c_pass"]) else "fail"
    rep[key] = res


# ---------------------------------------------------------------- H3
def mark_H3(df, rep):
    h = df[df["ins_A"].notna() & df["del"].notna()].copy()
    h["asym"] = h["ins_A"] - h["del"]
    out = {}
    for s, g in h.groupby("ss"):
        out[str(s)] = {"n": int(len(g)), "mean_asym": round(float(g["asym"].mean()), 3),
                       "ci95_by_domain": boot_by_domain(g.reset_index(drop=True),
                                                        lambda x: float(x["asym"].mean()), n=800)}
    rep["H3_exploratory"] = {"note": "exploratory, declared in advance; positive = insertion worse",
                             "by_ss": out}


# ---------------------------------------------------------------- permutation
def permutation(df, rep, seed=1):
    rng = np.random.default_rng(seed)
    p = df.copy()
    for col in ("ins_A", "del"):
        p[col] = p.groupby("domain")[col].transform(lambda s: rng.permutation(s.values))
    sub = {}
    mark_H1(p, sub, "ins_A", "H1_perm")
    mark_H2(p, sub, "ins_A", "H2_perm")
    rep["permutation_control"] = {
        "H1_difference_kcal": sub["H1_perm"]["difference_kcal"],
        "H1_ci95": sub["H1_perm"]["ci95_by_domain"],
        "H1_vanishes": bool(sub["H1_perm"]["ci95_by_domain"][0] <= 0 <= sub["H1_perm"]["ci95_by_domain"][1]),
        "H2_rho": sub["H2_perm"]["spearman_rho_asym_vs_dloop"],
        "H2_ci95": sub["H2_perm"]["ci95_by_domain"],
        "H2_vanishes": bool(sub["H2_perm"]["ci95_by_domain"][0] <= 0 <= sub["H2_perm"]["ci95_by_domain"][1])}


def main(rdir, out_path):
    df = load(rdir)
    rep = {"n_rows": int(len(df)), "n_domains": int(df["domain"].nunique()),
           "n_ins_A": int(df["ins_A"].notna().sum()), "n_ins_G": int(df["ins_G"].notna().sum()),
           "n_del": int(df["del"].notna().sum())}
    mark_B(df, rep)
    mark_H1(df, rep, "ins_A", "H1_ala")
    mark_H2(df, rep, "ins_A", "H2_ala")
    permutation(df, rep)
    mark_H3(df, rep)
    mark_H1(df, rep, "ins_G", "H1_gly_replication")
    mark_H2(df, rep, "ins_G", "H2_gly_replication")
    print(json.dumps(rep, indent=2))
    json.dump(rep, open(out_path, "w"), indent=2)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
