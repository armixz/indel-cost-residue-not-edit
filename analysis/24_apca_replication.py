"""Marks R1-R3 of PREREG_bulge.md: does the strand asymmetry replicate in an orthogonal assay?

Topolska et al. yeast abundancePCA (Zenodo 10830272), 9 domains, 1,450 deletions and 1,420
insertions, with the authors' own DSSP labels and rSASA. scaled_fitness is their scaled abundance
(higher = more abundant). cost = -scaled_fitness, so larger cost = worse, matching the Tsuboyama
convention used elsewhere in this project.

A = cost(insertion) - cost(deletion) at the same position. Negative = the insertion is cheaper.

Usage: indel_env python 24_apca_replication.py <apca_csv> <out_json>
"""
import json
import sys

import numpy as np
import pandas as pd

RNG = np.random.default_rng(0)
STRAND = {"Strand", "Bridge"}
HELIX = {"AlphaHelix", "310Helix", "PiHelix"}


def agg(vals, nboot=20000):
    v = np.asarray([x for x in vals if x is not None and np.isfinite(x)], dtype=float)
    if len(v) < 3:
        return {"n_domains": int(len(v)), "mean": None, "ci95": [None, None]}
    idx = RNG.integers(0, len(v), size=(nboot, len(v)))
    bs = v[idx].mean(axis=1)
    return {"n_domains": int(len(v)), "mean": round(float(v.mean()), 3),
            "values": [round(float(x), 3) for x in v],
            "ci95": [round(float(np.quantile(bs, 0.025)), 3),
                     round(float(np.quantile(bs, 0.975)), 3)]}


def build(path):
    d = pd.read_csv(path)
    d = d[d["mut_type"].isin(["insertions", "deletions"])].copy()
    d["cost"] = -d["scaled_fitness"]
    d["cls"] = np.where(d["Structure_name"].isin(STRAND), "E",
                        np.where(d["Structure_name"].isin(HELIX), "H", "L"))
    ins = (d[d["mut_type"] == "insertions"]
           .groupby(["domain", "Pos"])
           .agg(ins=("cost", "mean"), n_ins=("cost", "size")).reset_index())
    dele = (d[d["mut_type"] == "deletions"]
            .groupby(["domain", "Pos"])
            .agg(dele=("cost", "mean")).reset_index())
    meta = (d.groupby(["domain", "Pos"])
            .agg(cls=("cls", "first"), struct=("Structure_name", "first"),
                 rsasa=("rSASA", "first"), wt=("WT_aa", "first")).reset_index())
    m = ins.merge(dele, on=["domain", "Pos"]).merge(meta, on=["domain", "Pos"])
    m["asym"] = m["ins"] - m["dele"]
    return d, m


def per_domain_strand_minus_other(m):
    out = {}
    for dom, g in m.groupby("domain"):
        e, o = g[g["cls"] == "E"]["asym"], g[g["cls"] != "E"]["asym"]
        if len(e) >= 3 and len(o) >= 3:
            out[dom] = float(e.mean() - o.mean())
    return out


def main(apca_csv, out_path):
    raw, m = build(apca_csv)
    rep = {"n_positions_matched": int(len(m)), "n_domains": int(m["domain"].nunique()),
           "class_counts": {k: int(v) for k, v in m["cls"].value_counts().items()},
           "insertions_per_position_median": float(m["n_ins"].median())}

    rep["asymmetry_by_class"] = {str(c): agg([float(x["asym"].mean())
                                              for _, x in g.groupby("domain")])
                                 for c, g in m.groupby("cls")}

    # R1
    per = per_domain_strand_minus_other(m)
    neg = sum(1 for v in per.values() if v < 0)
    a = agg(list(per.values()))
    rep["R1_strand_minus_other"] = {
        "per_domain": {k: round(v, 3) for k, v in per.items()},
        "n_domains": len(per), "n_negative": neg, "summary": a,
        "mark_R1_pass": bool(neg >= 7 and a["ci95"][1] is not None and a["ci95"][1] < 0)}

    # R2 exposed only
    med = m["rsasa"].median()
    exp = m[m["rsasa"] > med]
    pere = per_domain_strand_minus_other(exp)
    ae = agg(list(pere.values()))
    rep["R2_exposed_only"] = {
        "rsasa_median": float(med), "n_positions": int(len(exp)),
        "per_domain": {k: round(v, 3) for k, v in pere.items()},
        "n_negative": int(sum(1 for v in pere.values() if v < 0)),
        "summary": ae,
        "by_class": {str(c): agg([float(x["asym"].mean()) for _, x in g.groupby("domain")])
                     for c, g in exp.groupby("cls")},
        "mark_R2_pass": bool(ae["ci95"][1] is not None and ae["ci95"][1] < 0)}

    # R3 identity control: inserted residue equals the wild-type residue at that position
    ins_rows = raw[raw["mut_type"] == "insertions"].copy()
    key = "insID" if ins_rows["insID"].notna().any() else "Mut"
    matched = ins_rows[ins_rows[key] == ins_rows["WT_aa"]]
    if len(matched) < 100:
        rep["R3_identity"] = {"n": int(len(matched)), "status": "UNTESTED - fewer than 100"}
        rep["mark_R3_pass"] = None
    else:
        mi = (matched.groupby(["domain", "Pos"])["cost"].mean().rename("ins_m").reset_index())
        mm = mi.merge(m[["domain", "Pos", "dele", "cls", "rsasa"]], on=["domain", "Pos"])
        mm["asym"] = mm["ins_m"] - mm["dele"]
        perm = {}
        for dom, g in mm.groupby("domain"):
            e, o = g[g["cls"] == "E"]["asym"], g[g["cls"] != "E"]["asym"]
            if len(e) >= 3 and len(o) >= 3:
                perm[dom] = float(e.mean() - o.mean())
        am = agg(list(perm.values()))
        rep["R3_identity"] = {"n_positions": int(len(mm)), "key_used": key,
                              "per_domain": {k: round(v, 3) for k, v in perm.items()},
                              "n_negative": int(sum(1 for v in perm.values() if v < 0)),
                              "summary": am,
                              "by_class": {str(c): agg([float(x["asym"].mean())
                                                        for _, x in g.groupby("domain")])
                                           for c, g in mm.groupby("cls")}}
        rep["mark_R3_pass"] = bool(am["mean"] is not None and am["mean"] < 0)

    # measurement error check (kill condition 3)
    err = raw.groupby("mut_type")["scaled_sigma"].describe()[["count", "mean", "50%"]]
    rep["error_check"] = json.loads(err.to_json(orient="index"))

    json.dump(rep, open(out_path, "w"), indent=2)
    print(json.dumps(rep, indent=2))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
