"""Marks D1-D4 of PREREG_decomposition.md.

Is a deletion an insertion plus the loss of a side chain?
  d = deletion cost, i = Ala insertion cost, g = cost of substituting wild-type -> Gly,
  all at the same position in the same protein.

Usage: indel_env python 37_decomp.py <lib3_csv> <results_dir> <afdir> <out_json>
"""
import csv
import json
import os
import re
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
exec(open(os.path.join(HERE, "34_scale.py")).read().split("def main")[0])  # noqa: S102

SUB = re.compile(r"^(.*?)_([A-Z])(\d+)([A-Z])$")
IND = re.compile(r"^(.*?)_(ins|del)([A-Z])(\d+)$")


def read_all(lib3):
    """Single indels and single substitutions from Lib3, with their dG and CI."""
    csv.field_size_limit(10 ** 7)
    wt, ind, sub = {}, [], []
    with open(lib3, newline="") as fh:
        for r in csv.DictReader(fh):
            nm = r["name"]
            try:
                dg = float(r["deltaG"])
                ci = float(r["deltaG_95CI"])
            except (ValueError, KeyError):
                continue
            tail = nm.split("_")[-1]
            if tail.startswith("wt") or tail.startswith("syn"):
                wt.setdefault(nm.rsplit("_", 1)[0], []).append(dg)
                continue
            m = IND.match(nm)
            if m:
                ind.append({"domain": m.group(1), "kind": m.group(2), "res": m.group(3),
                            "pos": int(m.group(4)), "dG": dg, "ci": ci})
                continue
            m = SUB.match(nm)
            if m:
                sub.append({"domain": m.group(1), "wtaa": m.group(2), "pos": int(m.group(3)),
                            "mut": m.group(4), "dG": dg, "ci": ci})
    return ({k: float(np.mean(v)) for k, v in wt.items()},
            pd.DataFrame(ind), pd.DataFrame(sub))


def build(lib3, afdir, prefix_filter=None):
    wt, ind, sub = read_all(lib3)
    if prefix_filter is not None:
        ind = ind[ind["domain"].str.startswith(prefix_filter)]
        sub = sub[sub["domain"].str.startswith(prefix_filter)]
    else:
        ind = ind[ind["domain"].str.endswith(".pdb")]
        sub = sub[sub["domain"].str.endswith(".pdb")]
    ind, sub = ind[ind["ci"] <= 1.0].copy(), sub[sub["ci"] <= 1.0].copy()
    for d in (ind, sub):
        d["wt"] = d["domain"].map(wt)
    ind, sub = ind[ind["wt"].notna()], sub[sub["wt"].notna()]
    ind = ind.assign(cost=ind["wt"] - ind["dG"])
    sub = sub.assign(cost=sub["wt"] - sub["dG"])

    i = (ind[(ind["kind"] == "ins") & (ind["res"] == "A")]
         .groupby(["domain", "pos"])["cost"].mean().rename("i").reset_index())
    d = (ind[ind["kind"] == "del"]
         .groupby(["domain", "pos"])["cost"].mean().rename("d").reset_index())
    g = (sub[sub["mut"] == "G"]
         .groupby(["domain", "pos"])["cost"].mean().rename("g").reset_index())
    wtaa = sub.groupby(["domain", "pos"])["wtaa"].first().rename("wtaa").reset_index()
    m = i.merge(d, on=["domain", "pos"]).merge(g, on=["domain", "pos"]).merge(
        wtaa, on=["domain", "pos"])

    # structure
    ssrows = []
    for dom in sorted(m["domain"].unique()):
        p = os.path.join(afdir, dom.replace(".pdb", "") + ".pdb")
        if not os.path.exists(p):
            continue
        res = d17.parse_backbone(p)
        s, keys = d17.assign(res)
        if not s:
            continue
        el = d17.elements(s, keys)
        for k in keys:
            ssrows.append({"domain": dom, "pos": k, "ss": el[k][0], "aa_model": res[k]["aa"]})
    S = pd.DataFrame(ssrows)
    if len(S):
        m = m.merge(S, on=["domain", "pos"], how="left")
    return m


def fit(df, xcols, ycol="d", nmin=12):
    a_s, b_s, rmse = [], [], []
    for _, grp in df.groupby("domain"):
        q = grp.dropna(subset=xcols + [ycol])
        if len(q) < nmin:
            continue
        X = np.column_stack([np.ones(len(q))] + [q[c].to_numpy() for c in xcols])
        if np.linalg.matrix_rank(X) < X.shape[1]:
            continue
        beta, *_ = np.linalg.lstsq(X, q[ycol].to_numpy(), rcond=None)
        pred = X @ beta
        a_s.append(float(beta[0]))
        b_s.append([float(x) for x in beta[1:]])
        rmse.append(float(np.sqrt(np.mean((q[ycol].to_numpy() - pred) ** 2))))
    return a_s, b_s, rmse


def main(lib3, rdir, afdir, out_path):
    global d17
    import importlib.util
    sp = importlib.util.spec_from_file_location("d17", os.path.join(HERE, "17_dssp.py"))
    d17 = importlib.util.module_from_spec(sp)
    sp.loader.exec_module(d17)

    rep = {}
    for label, pf in (("natural", None), ("designed", "HEEH")):
        m = build(lib3, afdir, pf)
        if len(m) == 0:
            rep[label] = {"status": "no data"}
            continue
        m["sum"] = m["i"] + m["g"]
        r = {"n_positions": int(len(m)), "n_domains": int(m["domain"].nunique()),
             "mean": {"i": round(float(m["i"].mean()), 3), "d": round(float(m["d"].mean()), 3),
                      "g": round(float(m["g"].mean()), 3)}}
        if len(m) < 3000 or m["domain"].nunique() < 100:
            r["power_note"] = "UNDERPOWERED by the pre-registered threshold"

        # D1
        a1, b1, e1 = fit(m, ["i"])
        a2, b2, e2 = fit(m, ["sum"])
        a3, b3, e3 = fit(m, ["i", "g"])
        r["d_on_i"] = {"intercept": agg(a1), "slope": agg([x[0] for x in b1]),
                       "rmse": round(float(np.mean(e1)), 3)}
        r["d_on_i_plus_g"] = {"intercept": agg(a2), "slope": agg([x[0] for x in b2]),
                              "rmse": round(float(np.mean(e2)), 3)}
        r["d_on_i_and_g_separately"] = {"intercept": agg(a3),
                                        "coef_i": agg([x[0] for x in b3]),
                                        "coef_g": agg([x[1] for x in b3]),
                                        "rmse": round(float(np.mean(e3)), 3)}
        drop = 1 - np.mean(e2) / np.mean(e1)
        sl = r["d_on_i_plus_g"]["slope"]
        r["mark_D1_pass"] = bool(sl["ci95"][0] is not None and sl["ci95"][0] <= 1.0 <= sl["ci95"][1]
                                 and drop >= 0.15)
        r["rmse_drop_fraction"] = round(float(drop), 3)

        # D2
        ic1, ic2 = r["d_on_i"]["intercept"], r["d_on_i_plus_g"]["intercept"]
        r["mark_D2_pass"] = bool(ic2["mean"] is not None and ic1["mean"] is not None
                                 and abs(ic2["mean"]) < abs(ic1["mean"])
                                 and not (ic2["ci95"][0] <= ic1["mean"] <= ic2["ci95"][1]))

        # D3 glycine positions
        gly = m[m["wtaa"] == "G"]
        allv = agg([float((x["d"] - x["i"]).mean()) for _, x in m.groupby("domain") if len(x) >= 3])
        glyv = agg([float((x["d"] - x["i"]).mean()) for _, x in gly.groupby("domain")
                    if len(x) >= 2])
        r["D3_glycine_positions"] = {"n": int(len(gly)), "offset_at_gly": glyv,
                                     "offset_overall": allv}
        r["mark_D3_pass"] = bool(glyv["mean"] is not None and abs(glyv["mean"]) < 0.2)

        # D4 strand excess before and after removing g
        if "ss" in m:
            m["asym"] = m["i"] - m["d"]
            m["asym_adj"] = m["i"] - (m["d"] - m["g"])
            def el_pair(col):
                return agg([float(x[x["ss"] == "E"][col].mean() - x[x["ss"] == "L"][col].mean())
                            for _, x in m.groupby("domain")
                            if (x["ss"] == "E").sum() >= 2 and (x["ss"] == "L").sum() >= 2])
            r["D4_strand_minus_loop"] = {"raw": el_pair("asym"), "side_chain_removed": el_pair("asym_adj")}
            adj = r["D4_strand_minus_loop"]["side_chain_removed"]["mean"]
            r["mark_D4_pass"] = bool(adj is not None and adj > -0.162)

        # control 3: permute g
        rng = np.random.default_rng(5)
        p = m.copy()
        p["g"] = p.groupby("domain")["g"].transform(lambda s: rng.permutation(s.values))
        p["sum"] = p["i"] + p["g"]
        _, _, ep = fit(p, ["sum"])
        r["permuted_g_rmse"] = round(float(np.mean(ep)), 3)
        rep[label] = r

    json.dump(rep, open(out_path, "w"), indent=2)
    print(json.dumps(rep, indent=2))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
