"""Adjudicating H-strand: why does the yeast abundance assay not reproduce it?

Three questions, in order of how decisive they are.

1. PERMUTATION STABILITY. The single permutation in 30_strand_controls.py gave +0.109
   [0.002, 0.225] for insertions, which does not cleanly contain zero, though its sign is opposite
   to the effect. Run 200 independent permutations and report the null distribution. If the null
   mean is reliably positive, the observed -0.478 is conservative; if it wanders, the earlier CI
   was a fluke of one draw.

2. SHARED DOMAINS. The nine abundancePCA domains may also exist in Tsuboyama Lib3. If a domain is
   in both, the two assays can be compared on the SAME protein, which separates "the effect is
   specific to thermodynamic stability" from "the effect is an artefact of one analysis".

3. DYNAMIC RANGE. An abundance assay saturates: a variant that is fully degraded cannot get worse.
   If strand indels sit at the floor more often than helix indels, the abundance assay cannot
   resolve the difference even if it exists. Measured as the fraction of variants at or below the
   5th percentile of the domain's own distribution, and by repeating the strand-helix comparison
   after dropping floored variants.

Usage: indel_env python 31_adjudicate.py <results_dir> <apca_csv> <out_json>
"""
import json
import sys

import numpy as np
import pandas as pd

CI_MAX = 1.0
RNG = np.random.default_rng(0)


def agg(vals, nboot=20000):
    v = np.asarray([x for x in vals if x is not None and np.isfinite(x)], dtype=float)
    if len(v) < 3:
        return {"n_domains": int(len(v)), "mean": None, "ci95": [None, None]}
    idx = RNG.integers(0, len(v), size=(nboot, len(v)))
    bs = v[idx].mean(axis=1)
    return {"n_domains": int(len(v)), "mean": round(float(v.mean()), 3),
            "ci95": [round(float(np.quantile(bs, 0.025)), 3),
                     round(float(np.quantile(bs, 0.975)), 3)]}


def paired_EH(df, col, nmin=2, sscol="ss", dom="domain"):
    out = []
    for _, g in df.groupby(dom):
        e, h = g[g[sscol] == "E"][col], g[g[sscol] == "H"][col]
        if len(e) >= nmin and len(h) >= nmin:
            out.append(float(e.mean() - h.mean()))
    return out


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
    return df[df["ins_A"].notna() & df["del"].notna()].copy()


def main(rdir, apca_csv, out_path):
    df = load(rdir)
    rep = {}

    # 1. permutation stability
    nulls_i, nulls_d = [], []
    for seed in range(200):
        rng = np.random.default_rng(1000 + seed)
        p = df.copy()
        for c in ("ins_A", "del"):
            p[c] = p.groupby("domain")[c].transform(lambda s: rng.permutation(s.values))
        vi = paired_EH(p, "ins_A")
        vd = paired_EH(p, "del")
        if vi:
            nulls_i.append(float(np.mean(vi)))
        if vd:
            nulls_d.append(float(np.mean(vd)))
    obs_i = float(np.mean(paired_EH(df, "ins_A")))
    obs_d = float(np.mean(paired_EH(df, "del")))
    rep["permutation_200"] = {
        "observed_insertion": round(obs_i, 3), "observed_deletion": round(obs_d, 3),
        "null_insertion": {"mean": round(float(np.mean(nulls_i)), 4),
                           "sd": round(float(np.std(nulls_i)), 4),
                           "q025": round(float(np.quantile(nulls_i, 0.025)), 3),
                           "q975": round(float(np.quantile(nulls_i, 0.975)), 3)},
        "null_deletion": {"mean": round(float(np.mean(nulls_d)), 4),
                          "sd": round(float(np.std(nulls_d)), 4),
                          "q025": round(float(np.quantile(nulls_d, 0.025)), 3),
                          "q975": round(float(np.quantile(nulls_d, 0.975)), 3)},
        "p_insertion": float(np.mean([x <= obs_i for x in nulls_i])),
        "p_deletion": float(np.mean([x <= obs_d for x in nulls_d]))}

    # 2. shared domains
    d = pd.read_csv(apca_csv)
    seqs = d.groupby("domain")["aa_seq"].first().to_dict()
    tsub_seq = (df.sort_values(["domain", "struct_pos"])
                .groupby("domain")["aa"].apply(lambda s: "".join(s)).to_dict())
    matches = {}
    for name, s in seqs.items():
        s = str(s)
        best = None
        for tname, tseq in tsub_seq.items():
            if len(tseq) < 15:
                continue
            if tseq in s or s in tseq:
                best = (tname, "containment", len(tseq))
                break
            # 20-mer anchor
            if any(tseq[i:i + 20] in s for i in range(0, max(1, len(tseq) - 20), 5)):
                best = (tname, "anchor20", len(tseq))
                break
        matches[name] = best
    rep["shared_domains"] = {k: (list(v) if v else None) for k, v in matches.items()}

    # 3. dynamic range in the abundance assay
    a = d[d["mut_type"].isin(["insertions", "deletions"])].copy()
    a["ss"] = np.where(a["Structure_name"].isin(["Strand", "Bridge"]), "E",
                       np.where(a["Structure_name"].isin(["AlphaHelix", "310Helix"]), "H", "L"))
    a["cost"] = -a["scaled_fitness"]
    # the floor: per domain, the 5th percentile of ALL variants of that domain
    floors = a.groupby("domain")["cost"].quantile(0.95).to_dict()   # cost high = dead
    a["floored"] = a.apply(lambda r: r["cost"] >= floors[r["domain"]] * 0.95, axis=1)
    rep["dynamic_range"] = {
        "floored_fraction_by_class": {str(c): round(float(g["floored"].mean()), 3)
                                      for c, g in a.groupby("ss")},
        "cost_sd_by_class": {str(c): round(float(g["cost"].std()), 3)
                             for c, g in a.groupby("ss")},
        "n_by_class": {str(c): int(len(g)) for c, g in a.groupby("ss")}}
    pos = (a[~a["floored"]].groupby(["domain", "Pos", "mut_type"])
           .agg(cost=("cost", "mean"), ss=("ss", "first")).reset_index())
    for mt, key in (("insertions", "insertion"), ("deletions", "deletion")):
        v = paired_EH(pos[pos["mut_type"] == mt], "cost")
        rep["dynamic_range"][f"EH_unfloored_{key}"] = agg(v)
        rep["dynamic_range"][f"EH_unfloored_{key}_n_negative"] = int(sum(1 for x in v if x < 0))

    # also: does the Tsuboyama effect hold in BETA-RICH domains only, which is what the aPCA set is?
    comp = df.groupby("domain")["ss"].value_counts(normalize=True).unstack(fill_value=0.0)
    for c in ("H", "E", "L"):
        if c not in comp:
            comp[c] = 0.0
    beta_rich = set(comp[comp["E"] >= comp["H"]].index)
    helix_rich = set(comp[comp["H"] > comp["E"]].index)
    rep["by_fold_bias"] = {
        "beta_biased_domains": {"n": len(beta_rich),
                                "insertion": agg(paired_EH(df[df["domain"].isin(beta_rich)], "ins_A")),
                                "deletion": agg(paired_EH(df[df["domain"].isin(beta_rich)], "del"))},
        "helix_biased_domains": {"n": len(helix_rich),
                                 "insertion": agg(paired_EH(df[df["domain"].isin(helix_rich)], "ins_A")),
                                 "deletion": agg(paired_EH(df[df["domain"].isin(helix_rich)], "del"))}}

    json.dump(rep, open(out_path, "w"), indent=2)
    print(json.dumps(rep, indent=2))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])
