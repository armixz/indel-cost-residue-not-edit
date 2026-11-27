"""Marks RS3/RS4 of PREREG_strand.md: H-strand in designed HEEH_KT miniproteins.

These designs carry single-residue indels in Lib3, have AlphaFold models here, and are built with
a helix-strand-strand-helix topology, so each one supplies a within-protein helix-versus-strand
contrast. They are not natural proteins and share no evolutionary history with the domains that
produced the observation.

Numbering is verified from the deletion variants, whose name carries the deleted wild-type
residue: the model residue at that position must match. The set is dropped if the match rate is
below 0.95.

Usage: indel_env python 33_designed.py <lib3_csv> <results_dir> <afdir> <out_json>
"""
import csv
import glob
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
        return {"n_domains": int(len(v)), "mean": None, "ci95": [None, None]}
    idx = RNG.integers(0, len(v), size=(nboot, len(v)))
    bs = v[idx].mean(axis=1)
    return {"n_domains": int(len(v)), "mean": round(float(v.mean()), 3),
            "ci95": [round(float(np.quantile(bs, 0.025)), 3),
                     round(float(np.quantile(bs, 0.975)), 3)]}


def paired_EH(df, col, nmin=2):
    out = []
    for _, g in df.groupby("domain"):
        e, h = g[g["ss"] == "E"][col], g[g["ss"] == "H"][col]
        if len(e) >= nmin and len(h) >= nmin:
            out.append(float(e.mean() - h.mean()))
    return out


def main(lib3, rdir, afdir, out_path):
    # wild-type dG and variants for designed proteins
    csv.field_size_limit(10 ** 7)
    wt, rows = {}, []
    with open(lib3, newline="") as fh:
        for r in csv.DictReader(fh):
            nm = r["name"]
            if not nm.startswith("HEEH"):
                continue
            if nm.split("_")[-1].startswith("wt"):
                # several wt rows per design (wth/wtm/wty/wte); keep their mean
                try:
                    wt.setdefault(nm.rsplit("_", 1)[0], []).append(float(r["deltaG"]))
                except ValueError:
                    pass
                continue
            m = PAT.match(nm)
            if not m:
                continue
            try:
                dg = float(r["deltaG"])
                ci = float(r["deltaG_95CI"])
            except (ValueError, KeyError):
                continue
            rows.append({"domain": m.group(1), "kind": m.group(2), "res": m.group(3),
                         "pos": int(m.group(4)), "dG": dg, "ci": ci})
    v = pd.DataFrame(rows)

    # structures
    ssrows, checks = [], []
    for dom in sorted(v["domain"].unique()):
        p = os.path.join(afdir, dom + ".pdb")
        if not os.path.exists(p):
            continue
        res = d17.parse_backbone(p)
        ss, keys = d17.assign(res)
        if not ss:
            continue
        el = d17.elements(ss, keys)
        a, b = min(keys), max(keys)
        for i in keys:
            ssrows.append({"domain": dom, "pos": i, "aa": res[i]["aa"], "ss": el[i][0],
                           "idx_end": el[i][4] if el[i][0] != "L" else -1,
                           "plddt": float(np.mean(res[i].get("b", [np.nan]))
                                          if isinstance(res[i].get("b"), list) else np.nan)})
        # numbering check with deletions
        d = v[(v["domain"] == dom) & (v["kind"] == "del")]
        for _, r in d.iterrows():
            if a <= r["pos"] <= b:
                checks.append(res[r["pos"]]["aa"] == r["res"])
    S = pd.DataFrame(ssrows)
    match = float(np.mean(checks)) if checks else 0.0

    rep = {"designs_with_variants": int(v["domain"].nunique()),
           "designs_with_models": int(S["domain"].nunique()) if len(S) else 0,
           "numbering_match_rate": round(match, 4),
           "n_checks": len(checks)}
    if match < 0.95:
        rep["status"] = "DROPPED - numbering match below 0.95"
        json.dump(rep, open(out_path, "w"), indent=2)
        print(json.dumps(rep, indent=2))
        return

    # wild-type reference: median dG of that design's own synonymous/wt rows, else fall back
    # to the design's median variant dG is NOT valid, so require a wt row
    v = v[v["ci"] <= 1.0]
    wtm = {k: float(sum(x) / len(x)) for k, x in wt.items()}
    v["wt"] = v["domain"].map(wtm)
    miss = int(v["wt"].isna().sum())
    v = v[v["wt"].notna()].copy()
    v["cost"] = v["wt"] - v["dG"]
    rep["variants_used"] = int(len(v))
    rep["variants_without_wt_reference"] = miss

    ins = v[(v["kind"] == "ins")].groupby(["domain", "pos"])["cost"].mean().rename("ins").reset_index()
    dele = v[(v["kind"] == "del")].groupby(["domain", "pos"])["cost"].mean().rename("del").reset_index()
    m = S.merge(ins, on=["domain", "pos"], how="left").merge(dele, on=["domain", "pos"], how="left")

    rep["ss_fractions"] = {k: round(float(x), 3) for k, x in
                           m["ss"].value_counts(normalize=True).items()}
    rep["positions_with_both"] = int((m["ins"].notna() & m["del"].notna()).sum())

    both = m[m["ins"].notna() & m["del"].notna()]
    vi = paired_EH(both, "ins")
    vd = paired_EH(both, "del")
    rep["RS3_insertion"] = {"summary": agg(vi), "n_negative": int(sum(1 for x in vi if x < 0)),
                            "n": len(vi)}
    rep["RS4_deletion"] = {"summary": agg(vd), "n_negative": int(sum(1 for x in vd if x < 0)),
                           "n": len(vd)}
    rep["mark_RS3_pass"] = bool(len(vi) >= 30 and rep["RS3_insertion"]["summary"]["ci95"][1] is not None
                                and rep["RS3_insertion"]["summary"]["ci95"][1] < 0)
    rep["mark_RS4_pass"] = bool(len(vd) >= 30 and rep["RS4_deletion"]["summary"]["ci95"][1] is not None
                                and rep["RS4_deletion"]["summary"]["ci95"][1] < 0)
    rep["class_means"] = {str(s): {"insertion": round(float(g["ins"].mean()), 3),
                                   "deletion": round(float(g["del"].mean()), 3),
                                   "n": int(len(g))}
                          for s, g in both.groupby("ss")}
    rep["asymmetry_by_ss"] = {str(s): agg([float((x["ins"] - x["del"]).mean())
                                           for _, x in g.groupby("domain")])
                              for s, g in both.groupby("ss")}

    # permutation
    rng = np.random.default_rng(3)
    p = both.copy()
    for c in ("ins", "del"):
        p[c] = p.groupby("domain")[c].transform(lambda s: rng.permutation(s.values))
    rep["permutation"] = {"insertion": agg(paired_EH(p, "ins")),
                          "deletion": agg(paired_EH(p, "del"))}

    json.dump(rep, open(out_path, "w"), indent=2)
    print(json.dumps(rep, indent=2))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
