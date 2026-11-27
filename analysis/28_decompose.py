"""Decompose the strand asymmetry, and ask whether strands behave differently in all-beta and
mixed alpha/beta domains.

27_within_domain.py left two things open.

(a) The asymmetry A = cost(insertion) - cost(deletion) is a difference. It could arise because
    insertions are unusually CHEAP in strands (what H-bulge predicts: the insertion has an escape
    route) or because deletions are unusually EXPENSIVE there (a different story: the ladder
    cannot close a gap). The two are distinguished by looking at the components against a
    within-domain baseline, which is what this script does.

(b) Within-domain E - L was -0.324 over 227 domains while E - H was -0.033 over the 84 domains
    that contain both. Those cannot both hold unless strands in all-beta domains differ from
    strands in mixed domains. Tested here directly.

Baseline: for each domain, the mean cost of loop positions. Costs are expressed relative to it, so
every number is a within-domain quantity and no between-domain variance enters.

Usage: indel_env python 28_decompose.py <results_dir> <out_json>
"""
import json
import sys

import numpy as np
import pandas as pd

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
    for key, kind, res in (("ins_A", "ins", "A"), ("ins_G", "ins", "G"), ("del", "del", None)):
        sub = ind[ind["kind"] == kind]
        if res is not None:
            sub = sub[sub["res"] == res]
        t = sub.groupby(["domain", "pos"])["cost"].mean().rename(key).reset_index()
        df = df.merge(t, on=["domain", "pos"], how="left")
    df = df[df["ins_A"].notna() & df["del"].notna()].copy()
    df["asym"] = df["ins_A"] - df["del"]
    return df


def rel_to_loop(df, col):
    """Per domain: mean of col in each class minus mean of col in that domain's loops."""
    out = {"E": [], "H": []}
    for _, g in df.groupby("domain"):
        base = g[g["ss"] == "L"][col]
        if len(base) < 2:
            continue
        b = float(base.mean())
        for c in ("E", "H"):
            x = g[g["ss"] == c][col]
            if len(x) >= 2:
                out[c].append(float(x.mean()) - b)
    return out


def main(rdir, out_path):
    df = load(rdir)
    rep = {"n_positions": int(len(df)), "n_domains": int(df["domain"].nunique())}

    # (a) components relative to the domain's own loops
    rep["relative_to_loop"] = {}
    for col, label in (("ins_A", "insertion_cost"), ("del", "deletion_cost"),
                       ("asym", "asymmetry")):
        r = rel_to_loop(df, col)
        rep["relative_to_loop"][label] = {c: agg(v) for c, v in r.items()}

    # absolute class means, for reference
    rep["absolute_class_means"] = {
        str(s): {"insertion": round(float(g["ins_A"].mean()), 3),
                 "deletion": round(float(g["del"].mean()), 3),
                 "n": int(len(g))}
        for s, g in df.groupby("ss")}

    # (b) all-beta vs mixed domains
    comp = df.groupby("domain")["ss"].value_counts(normalize=True).unstack(fill_value=0.0)
    for c in ("H", "E", "L"):
        if c not in comp:
            comp[c] = 0.0
    allbeta = set(comp[(comp["H"] == 0) & (comp["E"] > 0)].index)
    mixed = set(comp[(comp["H"] > 0) & (comp["E"] > 0)].index)
    allalpha = set(comp[(comp["E"] == 0) & (comp["H"] > 0)].index)

    def paired_EL(sub):
        out = []
        for _, g in sub.groupby("domain"):
            x, y = g[g["ss"] == "E"]["asym"], g[g["ss"] == "L"]["asym"]
            if len(x) >= 2 and len(y) >= 2:
                out.append(float(x.mean() - y.mean()))
        return agg(out)

    def paired_HL(sub):
        out = []
        for _, g in sub.groupby("domain"):
            x, y = g[g["ss"] == "H"]["asym"], g[g["ss"] == "L"]["asym"]
            if len(x) >= 2 and len(y) >= 2:
                out.append(float(x.mean() - y.mean()))
        return agg(out)

    rep["by_domain_class"] = {
        "counts": {"all_beta": len(allbeta), "mixed": len(mixed), "all_alpha": len(allalpha)},
        "E_minus_L_in_all_beta": paired_EL(df[df["domain"].isin(allbeta)]),
        "E_minus_L_in_mixed": paired_EL(df[df["domain"].isin(mixed)]),
        "H_minus_L_in_mixed": paired_HL(df[df["domain"].isin(mixed)]),
        "H_minus_L_in_all_alpha": paired_HL(df[df["domain"].isin(allalpha)])}

    # loop asymmetry itself, by domain class -- is the BASELINE moving?
    rep["loop_asymmetry_by_domain_class"] = {
        k: agg([float(g[g["ss"] == "L"]["asym"].mean())
                for _, g in df[df["domain"].isin(s)].groupby("domain")
                if (g["ss"] == "L").sum() >= 2])
        for k, s in (("all_beta", allbeta), ("mixed", mixed), ("all_alpha", allalpha))}

    # and the absolute insertion / deletion costs by domain class
    rep["absolute_by_domain_class"] = {
        k: {"insertion": round(float(df[df["domain"].isin(s)]["ins_A"].mean()), 3),
            "deletion": round(float(df[df["domain"].isin(s)]["del"].mean()), 3),
            "n_positions": int(df["domain"].isin(s).sum())}
        for k, s in (("all_beta", allbeta), ("mixed", mixed), ("all_alpha", allalpha))}

    json.dump(rep, open(out_path, "w"), indent=2)
    print(json.dumps(rep, indent=2))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
