"""Is the strand excess a within-protein structural effect, or a between-protein one?

Small domains tend to be all-alpha or all-beta, so the domains that contribute strand positions
are largely not the domains that contribute helix positions. A class difference computed by
averaging class means over domains is then a comparison between different proteins, not between
different parts of the same protein. 26_mechanism.py's within-domain paired strand-minus-helix
came out at -0.049 [-0.193, 0.096], which is why this script exists.

Three questions:
  1. paired within-domain class differences, with the number of domains that can supply the pair
  2. does a domain's overall beta content predict its overall asymmetry? (the between-protein
     alternative)
  3. a two-way decomposition: position-level asymmetry regressed on class with DOMAIN FIXED
     EFFECTS, which removes every between-domain difference by construction

Usage: indel_env python 27_within_domain.py <results_dir> <out_json>
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


def pair(df, a, b, nmin=2):
    out = []
    for _, g in df.groupby("domain"):
        x, y = g[g["ss"] == a]["asym"], g[g["ss"] == b]["asym"]
        if len(x) >= nmin and len(y) >= nmin:
            out.append(float(x.mean() - y.mean()))
    return agg(out)


def main(rdir, out_path):
    df = load(rdir)
    rep = {"n_positions": int(len(df)), "n_domains": int(df["domain"].nunique())}

    # how mixed are these domains?
    comp = df.groupby("domain")["ss"].value_counts(normalize=True).unstack(fill_value=0.0)
    for c in ("H", "E", "L"):
        if c not in comp:
            comp[c] = 0.0
    rep["domain_composition"] = {
        "n_domains": int(len(comp)),
        "pure_alpha_no_strand": int((comp["E"] == 0).sum()),
        "pure_beta_no_helix": int((comp["H"] == 0).sum()),
        "mixed_both": int(((comp["E"] > 0) & (comp["H"] > 0)).sum()),
        "median_E_fraction": round(float(comp["E"].median()), 3),
        "median_H_fraction": round(float(comp["H"].median()), 3)}

    # 1. paired within-domain
    rep["paired_within_domain"] = {}
    for nmin in (1, 2, 3):
        rep["paired_within_domain"][f"min_{nmin}_positions_each"] = {
            "E_minus_H": pair(df, "E", "H", nmin),
            "E_minus_L": pair(df, "E", "L", nmin),
            "H_minus_L": pair(df, "H", "L", nmin)}

    # unpaired, for contrast
    rep["unpaired_class_means"] = {
        str(s): agg([float(x["asym"].mean()) for _, x in g.groupby("domain")])
        for s, g in df.groupby("ss")}

    # 2. between-protein alternative: beta content vs domain asymmetry
    dom = df.groupby("domain").agg(asym=("asym", "mean"), n=("asym", "size")).reset_index()
    dom = dom.merge(comp.reset_index()[["domain", "E", "H"]], on="domain")
    dom = dom[dom["n"] >= 10]
    rep["between_protein"] = {
        "n_domains": int(len(dom)),
        "spearman_beta_fraction_vs_asymmetry":
            round(float(stats.spearmanr(dom["E"], dom["asym"])[0]), 3),
        "p": float(stats.spearmanr(dom["E"], dom["asym"])[1]),
        "spearman_helix_fraction_vs_asymmetry":
            round(float(stats.spearmanr(dom["H"], dom["asym"])[0]), 3),
        "mean_asym_beta_rich_E_gt_0.3": round(float(dom[dom["E"] > 0.3]["asym"].mean()), 3),
        "n_beta_rich": int((dom["E"] > 0.3).sum()),
        "mean_asym_helix_rich_H_gt_0.3": round(float(dom[dom["H"] > 0.3]["asym"].mean()), 3),
        "n_helix_rich": int((dom["H"] > 0.3).sum())}

    # 3. domain fixed effects
    try:
        import statsmodels.formula.api as smf
        d = df.copy()
        d["dom"] = d["domain"].astype("category")
        m = smf.ols("asym ~ C(ss) + C(dom)", data=d).fit(
            cov_type="cluster", cov_kwds={"groups": d["dom"].cat.codes})
        keep = {k: round(float(v), 4) for k, v in m.params.items() if k.startswith("C(ss)")}
        cis = {k: [round(float(a), 4), round(float(b), 4)]
               for k, (a, b) in m.conf_int().iterrows() if k.startswith("C(ss)")}
        rep["domain_fixed_effects"] = {
            "n": int(len(d)), "coef": keep, "ci95": cis,
            "note": "reference class is E (alphabetically first); a POSITIVE coefficient for H or "
                    "L means that class is less asymmetric than strand WITHIN the same domain"}
        # same, restricted to domains that actually contain both E and H
        mixed = set(comp[(comp["E"] > 0) & (comp["H"] > 0)].index)
        dm = d[d["domain"].isin(mixed)].copy()
        dm["dom"] = dm["domain"].astype("category")
        m2 = smf.ols("asym ~ C(ss) + C(dom)", data=dm).fit(
            cov_type="cluster", cov_kwds={"groups": dm["dom"].cat.codes})
        rep["domain_fixed_effects_mixed_only"] = {
            "n": int(len(dm)), "n_domains": int(dm["domain"].nunique()),
            "coef": {k: round(float(v), 4) for k, v in m2.params.items() if k.startswith("C(ss)")},
            "ci95": {k: [round(float(a), 4), round(float(b), 4)]
                     for k, (a, b) in m2.conf_int().iterrows() if k.startswith("C(ss)")}}
    except ImportError:
        rep["domain_fixed_effects"] = {"status": "statsmodels unavailable"}

    json.dump(rep, open(out_path, "w"), indent=2)
    print(json.dumps(rep, indent=2))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
