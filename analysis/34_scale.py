"""Is the designed-protein replication a refutation of H-strand, or is it out of range?

Designed HEEH_KT miniproteins are marginally stable. Mean indel cost there is 0.42-0.55 kcal/mol
against 1.71-2.19 in the natural domains, so the dynamic range is about four times smaller and an
effect of -0.478 would appear as roughly -0.11, which sits at the edge of the observed CI
[-0.111, 0.067]. This script decides between "refuted" and "underpowered" by putting both sets on
one scale.

  1. wild-type dG distributions of the two sets
  2. the strand-minus-helix difference expressed as a FRACTION of each protein's wild-type dG
  3. the natural domains restricted to the designs' stability range, and the designs restricted to
     their most stable members
  4. the same comparison for the asymmetry, which did replicate, as a positive control that the
     designed set has usable signal

Usage: indel_env python 34_scale.py <lib3_csv> <results_dir> <afdir> <out_json>
"""
import csv
import importlib.util
import json
import os
import re
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("d17", os.path.join(HERE, "17_dssp.py"))
d17 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(d17)

RNG = np.random.default_rng(0)
PAT = re.compile(r"^(.*?)_(ins|del)([A-Z])(\d+)$")


def agg(vals, nboot=20000):
    v = np.asarray([x for x in vals if x is not None and np.isfinite(x)], dtype=float)
    if len(v) < 5:
        return {"n": int(len(v)), "mean": None, "ci95": [None, None]}
    idx = RNG.integers(0, len(v), size=(nboot, len(v)))
    bs = v[idx].mean(axis=1)
    return {"n": int(len(v)), "mean": round(float(v.mean()), 3),
            "ci95": [round(float(np.quantile(bs, 0.025)), 3),
                     round(float(np.quantile(bs, 0.975)), 3)]}


def pairs(df, col, scale=None, nmin=2):
    out = []
    for d, g in df.groupby("domain"):
        e, h = g[g["ss"] == "E"][col], g[g["ss"] == "H"][col]
        if len(e) >= nmin and len(h) >= nmin:
            v = float(e.mean() - h.mean())
            if scale is not None:
                s = scale.get(d)
                if s is None or s <= 0.3:
                    continue
                v = v / s
            out.append(v)
    return out


def designed(lib3, afdir):
    csv.field_size_limit(10 ** 7)
    wt, rows = {}, []
    with open(lib3, newline="") as fh:
        for r in csv.DictReader(fh):
            nm = r["name"]
            if not nm.startswith("HEEH"):
                continue
            if nm.split("_")[-1].startswith("wt"):
                try:
                    wt.setdefault(nm.rsplit("_", 1)[0], []).append(float(r["deltaG"]))
                except ValueError:
                    pass
                continue
            m = PAT.match(nm)
            if not m:
                continue
            try:
                rows.append({"domain": m.group(1), "kind": m.group(2), "pos": int(m.group(4)),
                             "dG": float(r["deltaG"]), "ci": float(r["deltaG_95CI"])})
            except (ValueError, KeyError):
                continue
    v = pd.DataFrame(rows)
    wtm = {k: float(np.mean(x)) for k, x in wt.items()}
    v = v[v["ci"] <= 1.0].copy()
    v["wt"] = v["domain"].map(wtm)
    v = v[v["wt"].notna()]
    v["cost"] = v["wt"] - v["dG"]
    ss = []
    for dom in sorted(v["domain"].unique()):
        p = os.path.join(afdir, dom + ".pdb")
        if not os.path.exists(p):
            continue
        res = d17.parse_backbone(p)
        s, keys = d17.assign(res)
        if not s:
            continue
        el = d17.elements(s, keys)
        for i in keys:
            ss.append({"domain": dom, "pos": i, "ss": el[i][0]})
    S = pd.DataFrame(ss)
    ins = v[v["kind"] == "ins"].groupby(["domain", "pos"])["cost"].mean().rename("ins").reset_index()
    de = v[v["kind"] == "del"].groupby(["domain", "pos"])["cost"].mean().rename("del").reset_index()
    m = S.merge(ins, on=["domain", "pos"]).merge(de, on=["domain", "pos"])
    m["asym"] = m["ins"] - m["del"]
    return m, wtm


def natural(rdir):
    ind = pd.read_csv(f"{rdir}/tsuboyama_indels.csv")
    feat = pd.read_csv(f"{rdir}/tsub_features.csv")
    ss = pd.read_csv(f"{rdir}/tsub_ss.csv")
    ss["domain"] = ss["domain"].str.replace(".pdb", "", regex=False)
    wt = json.load(open(f"{rdir}/tsub_wt_dG.json"))
    ind = ind[ind["dG_ci"] <= 1.0].copy()
    ind["wt"] = ind["domain"].map(wt)
    ind = ind[ind["wt"].notna()]
    ind["cost"] = ind["wt"] - ind["dG"]
    df = feat[["domain", "pos", "struct_pos"]].merge(
        ss[["domain", "struct_pos", "ss"]], on=["domain", "struct_pos"])
    for key, kind, res in (("ins", "ins", "A"), ("del", "del", None)):
        sub = ind[ind["kind"] == kind]
        if res is not None:
            sub = sub[sub["res"] == res]
        t = sub.groupby(["domain", "pos"])["cost"].mean().rename(key).reset_index()
        df = df.merge(t, on=["domain", "pos"], how="left")
    df = df[df["ins"].notna() & df["del"].notna()].copy()
    df["asym"] = df["ins"] - df["del"]
    return df, wt


def main(lib3, rdir, afdir, out_path):
    D, wtD = designed(lib3, afdir)
    N, wtN = natural(rdir)
    rep = {}

    dw = np.array([wtD[d] for d in D["domain"].unique() if d in wtD])
    nw = np.array([wtN[d] for d in N["domain"].unique() if d in wtN])
    rep["wild_type_dG"] = {
        "designed": {"n": len(dw), "median": round(float(np.median(dw)), 2),
                     "iqr": [round(float(np.quantile(dw, .25)), 2),
                             round(float(np.quantile(dw, .75)), 2)]},
        "natural": {"n": len(nw), "median": round(float(np.median(nw)), 2),
                    "iqr": [round(float(np.quantile(nw, .25)), 2),
                            round(float(np.quantile(nw, .75)), 2)]}}
    rep["mean_indel_cost"] = {
        "designed": {"insertion": round(float(D["ins"].mean()), 3),
                     "deletion": round(float(D["del"].mean()), 3)},
        "natural": {"insertion": round(float(N["ins"].mean()), 3),
                    "deletion": round(float(N["del"].mean()), 3)}}

    # 2. scaled by wild-type dG
    rep["EH_absolute"] = {"designed_ins": agg(pairs(D, "ins")), "natural_ins": agg(pairs(N, "ins")),
                          "designed_del": agg(pairs(D, "del")), "natural_del": agg(pairs(N, "del"))}
    rep["EH_scaled_by_wt_dG"] = {
        "designed_ins": agg(pairs(D, "ins", wtD)), "natural_ins": agg(pairs(N, "ins", wtN)),
        "designed_del": agg(pairs(D, "del", wtD)), "natural_del": agg(pairs(N, "del", wtN))}

    # 3. range matching
    lo, hi = float(np.quantile(dw, .05)), float(np.quantile(dw, .95))
    nmatch = N[N["domain"].map(lambda d: lo <= wtN.get(d, -9) <= hi)]
    rep["range_matched"] = {
        "designed_dG_range": [round(lo, 2), round(hi, 2)],
        "natural_domains_in_range": int(nmatch["domain"].nunique()),
        "natural_ins_in_that_range": agg(pairs(nmatch, "ins")),
        "natural_del_in_that_range": agg(pairs(nmatch, "del"))}
    top = sorted(wtD.items(), key=lambda kv: -kv[1])[:max(10, len(wtD) // 2)]
    topset = {k for k, _ in top}
    rep["most_stable_designs"] = {
        "n": len(topset),
        "insertion": agg(pairs(D[D["domain"].isin(topset)], "ins")),
        "deletion": agg(pairs(D[D["domain"].isin(topset)], "del"))}

    # 4. asymmetry as a positive control
    rep["asymmetry_positive_control"] = {
        "designed_by_ss": {str(s): agg([float(x["asym"].mean()) for _, x in g.groupby("domain")])
                           for s, g in D.groupby("ss")},
        "natural_by_ss": {str(s): agg([float(x["asym"].mean()) for _, x in g.groupby("domain")])
                          for s, g in N.groupby("ss")},
        "designed_EL_paired": agg([float(g[g["ss"] == "E"]["asym"].mean()
                                         - g[g["ss"] == "L"]["asym"].mean())
                                   for _, g in D.groupby("domain")
                                   if (g["ss"] == "E").sum() >= 2 and (g["ss"] == "L").sum() >= 2]),
        "natural_EL_paired": agg([float(g[g["ss"] == "E"]["asym"].mean()
                                        - g[g["ss"] == "L"]["asym"].mean())
                                  for _, g in N.groupby("domain")
                                  if (g["ss"] == "E").sum() >= 2 and (g["ss"] == "L").sum() >= 2])}

    json.dump(rep, open(out_path, "w"), indent=2)
    print(json.dumps(rep, indent=2))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
