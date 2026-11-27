"""Marks M1-M3 and the kill conditions of PREREG_bulge.md, on Tsuboyama Lib3.

M1 antiparallel vs parallel bridges
M2 edge (1 bridge partner) vs interior (2)
M3 strand middle vs strand end
Kill conditions: strand-length matching, permutation, measurement error, strand floor of 3.

All statistics are computed within domain and averaged over domains.

Usage: run.sh 26_mechanism.py <results_dir> <out_json>
"""
import json
import sys

import numpy as np
import pandas as pd

CI_MAX = 1.0
RNG = np.random.default_rng(0)


def agg(vals, nboot=10000):
    v = np.asarray([x for x in vals if x is not None and np.isfinite(x)], dtype=float)
    if len(v) < 10:
        return {"n_domains": int(len(v)), "mean": None, "ci95": [None, None]}
    idx = RNG.integers(0, len(v), size=(nboot, len(v)))
    bs = v[idx].mean(axis=1)
    return {"n_domains": int(len(v)), "mean": round(float(v.mean()), 3),
            "ci95": [round(float(np.quantile(bs, 0.025)), 3),
                     round(float(np.quantile(bs, 0.975)), 3)]}


def paired(df, mask_a, mask_b, col="asym"):
    """Within-domain mean(a) - mean(b), averaged over domains that have both."""
    out = []
    for _, g in df.groupby("domain"):
        a, b = g[mask_a(g)][col], g[mask_b(g)][col]
        if len(a) >= 1 and len(b) >= 1:
            out.append(float(a.mean() - b.mean()))
    return agg(out)


def load(rdir, strand_floor=2):
    ind = pd.read_csv(f"{rdir}/tsuboyama_indels.csv")
    feat = pd.read_csv(f"{rdir}/tsub_features.csv")
    ss = pd.read_csv(f"{rdir}/tsub_ss.csv")
    ss["domain"] = ss["domain"].str.replace(".pdb", "", regex=False)
    br = pd.read_csv(f"{rdir}/tsub_bridges.csv")
    wt = json.load(open(f"{rdir}/tsub_wt_dG.json"))

    ci = ind.groupby("kind")["dG_ci"].median().to_dict()
    ind = ind[ind["dG_ci"] <= CI_MAX].copy()
    ind["wt"] = ind["domain"].map(wt)
    ind = ind[ind["wt"].notna()]
    ind["cost"] = ind["wt"] - ind["dG"]

    df = feat[["domain", "pos", "struct_pos", "plddt", "rsa"]].merge(
        ss[["domain", "struct_pos", "aa", "ss", "idx_end", "elem_len"]],
        on=["domain", "struct_pos"], how="inner").merge(
        br, on=["domain", "struct_pos"], how="left")
    for key, kind, res in (("ins_A", "ins", "A"), ("del", "del", None)):
        sub = ind[ind["kind"] == kind]
        if res is not None:
            sub = sub[sub["res"] == res]
        t = sub.groupby(["domain", "pos"])["cost"].mean().rename(key).reset_index()
        df = df.merge(t, on=["domain", "pos"], how="left")
    df = df[df["ins_A"].notna() & df["del"].notna()].copy()
    df["asym"] = df["ins_A"] - df["del"]
    if strand_floor > 2:
        df.loc[(df["ss"] == "E") & (df["elem_len"] < strand_floor), "ss"] = "L"
    return df, ci


def main(rdir, out_path):
    df, ci_med = load(rdir)
    e = df[df["ss"] == "E"].copy()
    rep = {"strand_positions": int(len(e)), "strand_domains": int(e["domain"].nunique()),
           "median_dG_ci_by_kind": {k: round(float(v), 4) for k, v in ci_med.items()}}

    # M1 antiparallel vs parallel
    rep["M1_anti_minus_par"] = paired(e, lambda g: g["bridge_type"] == "anti",
                                      lambda g: g["bridge_type"] == "par")
    rep["M1_counts"] = {str(k): int(v) for k, v in e["bridge_type"].value_counts().items()}
    rep["M1_by_type"] = {str(t): agg([float(x["asym"].mean()) for _, x in g.groupby("domain")])
                         for t, g in e.groupby("bridge_type")}
    m1 = rep["M1_anti_minus_par"]
    rep["mark_M1_pass"] = bool(m1["mean"] is not None and m1["mean"] <= -0.15
                               and m1["ci95"][1] < 0)

    # M2 edge vs interior
    rep["M2_edge_minus_interior"] = paired(e, lambda g: g["n_partners"] == 1,
                                           lambda g: g["n_partners"] == 2)
    rep["M2_by_partners"] = {str(int(t)): agg([float(x["asym"].mean())
                                               for _, x in g.groupby("domain")])
                             for t, g in e[e["n_partners"].notna()].groupby("n_partners")}
    m2 = rep["M2_edge_minus_interior"]
    rep["mark_M2_pass"] = bool(m2["mean"] is not None and m2["mean"] <= -0.15
                               and m2["ci95"][1] < 0)

    # M3 strand middle vs strand end
    rep["M3_middle_minus_end"] = paired(e, lambda g: g["idx_end"] >= 2,
                                        lambda g: g["idx_end"] == 0)
    m3 = rep["M3_middle_minus_end"]
    rep["mark_M3_pass"] = bool(m3["mean"] is not None and m3["ci95"][1] < 0)

    # kill condition 1: element length matched
    len_match = {}
    for lo, hi in ((2, 4), (5, 6), (7, 20)):
        sub = df[(df["elem_len"] >= lo) & (df["elem_len"] <= hi)]
        len_match[f"len_{lo}_{hi}"] = {
            "n": int(len(sub)),
            "by_ss": {str(s): agg([float(x["asym"].mean()) for _, x in g.groupby("domain")])
                      for s, g in sub.groupby("ss")},
            "strand_minus_helix": paired(sub, lambda g: g["ss"] == "E",
                                         lambda g: g["ss"] == "H")}
    rep["kill1_length_matched"] = len_match

    # kill condition 2: permutation
    rng = np.random.default_rng(7)
    p = df.copy()
    p["asym"] = p.groupby("domain")["asym"].transform(lambda s: rng.permutation(s.values))
    pe = p[p["ss"] == "E"]
    rep["kill2_permutation"] = {
        "M1": paired(pe, lambda g: g["bridge_type"] == "anti", lambda g: g["bridge_type"] == "par"),
        "M2": paired(pe, lambda g: g["n_partners"] == 1, lambda g: g["n_partners"] == 2),
        "M3": paired(pe, lambda g: g["idx_end"] >= 2, lambda g: g["idx_end"] == 0),
        "strand_minus_helix": paired(p, lambda g: g["ss"] == "E", lambda g: g["ss"] == "H")}

    # kill condition 4: strand floor of 3
    df3, _ = load(rdir, strand_floor=3)
    rep["kill4_strand_floor_3"] = {
        "by_ss": {str(s): agg([float(x["asym"].mean()) for _, x in g.groupby("domain")])
                  for s, g in df3.groupby("ss")},
        "strand_minus_helix": paired(df3, lambda g: g["ss"] == "E", lambda g: g["ss"] == "H")}

    # headline for reference
    rep["headline_by_ss"] = {str(s): agg([float(x["asym"].mean()) for _, x in g.groupby("domain")])
                             for s, g in df.groupby("ss")}
    rep["headline_strand_minus_helix"] = paired(df, lambda g: g["ss"] == "E",
                                                lambda g: g["ss"] == "H")
    rep["H_bulge_verdict"] = ("supported" if (rep["mark_M1_pass"] and rep["mark_M2_pass"])
                              else "partial" if (rep["mark_M1_pass"] or rep["mark_M2_pass"])
                              else "REJECTED - both load-bearing marks failed")

    json.dump(rep, open(out_path, "w"), indent=2)
    print(json.dumps(rep, indent=2))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
