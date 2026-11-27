"""Controls C1-C5 and kill conditions of PREREG_strand.md, plus replication marks RS1/RS2.

H-strand: within a protein, a single-residue indel costs less in a beta strand than in an alpha
helix. Observed paired difference (84 mixed domains): insertions -0.478 [-0.665, -0.287],
deletions -0.445 [-0.669, -0.227].

Usage: indel_env python 30_strand_controls.py <results_dir> <apca_csv> <sec_struc_rds> <out_json>
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


def paired_EH(df, col, nmin=2, sscol="ss"):
    out = []
    for _, g in df.groupby("domain"):
        e, h = g[g[sscol] == "E"][col], g[g[sscol] == "H"][col]
        if len(e) >= nmin and len(h) >= nmin:
            out.append(float(e.mean() - h.mean()))
    return agg(out)


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
    return df[df["ins_A"].notna() & df["del"].notna()].copy()


def main(rdir, apca_csv, rds, out_path):
    df = load(rdir)
    rep = {"headline": {"insertion": paired_EH(df, "ins_A"),
                        "deletion": paired_EH(df, "del"),
                        "glycine_insertion": paired_EH(df[df["ins_G"].notna()], "ins_G")}}

    # C1 burial, C2 confidence
    for name, var in (("C1_burial_rsa", "rsa"), ("C2_confidence_plddt", "plddt")):
        s = df[df[var].notna()].copy()
        s["t"] = pd.qcut(s[var], 3, labels=["low", "mid", "high"], duplicates="drop")
        rep[name] = {str(t): {"insertion": paired_EH(g, "ins_A"), "deletion": paired_EH(g, "del"),
                              "n": int(len(g))}
                     for t, g in s.groupby("t", observed=True)}
        signs = [rep[name][t]["insertion"]["mean"] for t in rep[name]]
        rep[name + "_all_negative"] = bool(all(x is not None and x < 0 for x in signs))

    # C3 element position matched
    rep["C3_position_matched"] = {}
    for label, lo, hi in (("ends_idx0_1", 0, 1), ("middles_idx2plus", 2, 99)):
        s = df[(df["idx_end"] >= lo) & (df["idx_end"] <= hi) & (df["ss"].isin(["E", "H"]))]
        rep["C3_position_matched"][label] = {"insertion": paired_EH(s, "ins_A"),
                                             "deletion": paired_EH(s, "del"),
                                             "n": int(len(s))}
    rep["C3_pass"] = bool(all(v["insertion"]["mean"] is not None and v["insertion"]["mean"] < 0
                              for v in rep["C3_position_matched"].values()))

    # C4 residue composition
    per_aa = {}
    for a, g in df[df["ss"].isin(["E", "H"])].groupby("aa"):
        p = paired_EH(g, "ins_A", nmin=1)
        if p["mean"] is not None:
            per_aa[str(a)] = p
    neg = sum(1 for v in per_aa.values() if v["mean"] < 0)
    rep["C4_per_residue"] = {"per_aa": {k: v["mean"] for k, v in per_aa.items()},
                             "n_types": len(per_aa), "n_negative": neg,
                             "pass": bool(neg >= 0.7 * len(per_aa))}
    # and the same-residue restriction, pooled over residue types with domain pairing
    out = []
    for (_, a), g in df[df["ss"].isin(["E", "H"])].groupby(["domain", "aa"]):
        e, h = g[g["ss"] == "E"]["ins_A"], g[g["ss"] == "H"]["ins_A"]
        if len(e) >= 1 and len(h) >= 1:
            out.append(float(e.mean() - h.mean()))
    rep["C4_same_residue_pairs"] = agg(out)

    # C5 independent labels
    try:
        import pyreadr
        their = list(pyreadr.read_r(rds).values())[0]
        their = their.rename(columns={"Pos": "struct_pos", "pdb_name": "domain"})
        their["struct_pos"] = their["struct_pos"].astype(int)
        MAP = {"AlphaHelix": "H", "310Helix": "H", "PiHelix": "H", "Strand": "E", "Bridge": "E"}
        their["ss_them"] = their["simple_struc"].map(MAP).fillna("L")
        j = df.merge(their[["domain", "struct_pos", "ss_them"]], on=["domain", "struct_pos"])
        rep["C5_their_labels"] = {"n": int(len(j)),
                                  "insertion": paired_EH(j, "ins_A", sscol="ss_them"),
                                  "deletion": paired_EH(j, "del", sscol="ss_them"),
                                  "insertion_my_labels_same_rows": paired_EH(j, "ins_A")}
    except Exception as e:  # noqa: BLE001
        rep["C5_their_labels"] = {"error": str(e)}

    # kill condition 1: permutation
    rng = np.random.default_rng(11)
    p = df.copy()
    for c in ("ins_A", "del"):
        p[c] = p.groupby("domain")[c].transform(lambda s: rng.permutation(s.values))
    rep["kill_permutation"] = {"insertion": paired_EH(p, "ins_A"), "deletion": paired_EH(p, "del")}

    # ---- replication in aPCA
    d = pd.read_csv(apca_csv)
    d = d[d["mut_type"].isin(["insertions", "deletions"])].copy()
    d["cost"] = -d["scaled_fitness"]
    d["ss"] = np.where(d["Structure_name"].isin(["Strand", "Bridge"]), "E",
                       np.where(d["Structure_name"].isin(["AlphaHelix", "310Helix"]), "H", "L"))
    d = d.rename(columns={"Pos": "struct_pos"})
    pos = (d.groupby(["domain", "struct_pos", "mut_type"])
           .agg(cost=("cost", "mean"), ss=("ss", "first"), rsasa=("rSASA", "first")).reset_index())
    rep["RS_replication"] = {}
    for mt, key in (("insertions", "RS1_insertion"), ("deletions", "RS2_deletion")):
        s = pos[pos["mut_type"] == mt]
        per = {}
        for dom, g in s.groupby("domain"):
            e, h = g[g["ss"] == "E"]["cost"], g[g["ss"] == "H"]["cost"]
            if len(e) >= 2 and len(h) >= 2:
                per[dom] = float(e.mean() - h.mean())
        a = agg(list(per.values()))
        rep["RS_replication"][key] = {
            "per_domain": {k: round(v, 3) for k, v in per.items()},
            "n_domains": len(per), "n_negative": int(sum(1 for v in per.values() if v < 0)),
            "summary": a,
            "mark_pass": bool(len(per) > 0 and sum(1 for v in per.values() if v < 0) >= 6
                              and a["ci95"][1] is not None and a["ci95"][1] < 0)}
        # class means for context
        rep["RS_replication"][key]["class_means"] = {
            str(c): round(float(g["cost"].mean()), 3) for c, g in s.groupby("ss")}

    json.dump(rep, open(out_path, "w"), indent=2)
    print(json.dumps(rep, indent=2))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
