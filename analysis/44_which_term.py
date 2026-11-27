"""Does the propensity correlation live in the helix term or the loop term?

The reported quantity is a contrast: mean deletion cost in helices minus mean deletion cost in
loops, per residue. A correlation with helix propensity could arise from either side. If deleting
glycine from a LOOP were unusually cheap, the contrast for glycine would be large for a reason
that says nothing about helices, and the propensity correlation would be an artefact of the
baseline rather than a statement about helix fit.

This script splits the contrast:
  1. helix cost alone, expressed relative to each domain's overall mean deletion cost
  2. loop cost alone, same normalisation
  3. strand cost alone, as a third comparison
  4. correlation of each with helix propensity
  5. a within-residue z-score version, to confirm the result does not depend on the normalisation

The prediction of the propensity account is that the correlation sits in the HELIX term:
good helix formers are expensive to remove from helices, while their cost in loops carries no
particular ordering by helix propensity.

Usage: indel_env python 44_which_term.py <lib3_csv> <afdir> <out_json>
"""
import csv
import importlib.util
import json
import os
import re
import sys
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
sp = importlib.util.spec_from_file_location("d17", os.path.join(HERE, "17_dssp.py"))
D17 = importlib.util.module_from_spec(sp)
sp.loader.exec_module(D17)
PAT = re.compile(r"^(?P<dom>.+?)_(?P<kind>ins|del)(?P<res>[A-Z])(?P<pos>\d+)$")
RNG = np.random.default_rng(0)
PACE = {"A": 0.00, "L": 0.21, "R": 0.21, "M": 0.24, "K": 0.26, "Q": 0.39, "E": 0.40,
        "I": 0.41, "W": 0.49, "S": 0.50, "Y": 0.53, "F": 0.54, "H": 0.61, "V": 0.61,
        "N": 0.65, "T": 0.66, "C": 0.68, "D": 0.69, "G": 1.00}


def load(lib3, afdir):
    csv.field_size_limit(10 ** 7)
    wt, rows = {}, []
    with open(lib3, newline="") as fh:
        for r in csv.DictReader(fh):
            nm = r["name"]
            if nm.rsplit("_", 1)[-1].startswith("wt"):
                wt.setdefault(nm.rsplit("_", 1)[0], []).append(float(r["deltaG"]))
                continue
            m = PAT.match(nm)
            if not m or m["kind"] != "del":
                continue
            try:
                rows.append({"domain": m["dom"], "aa": m["res"], "pos": int(m["pos"]),
                             "dG": float(r["deltaG"]), "ci": float(r["deltaG_95CI"])})
            except (ValueError, KeyError):
                continue
    wtm = {k: float(np.mean(v)) for k, v in wt.items()}
    v = pd.DataFrame(rows)
    v = v[v["ci"] <= 1.0].copy()
    v["wtdG"] = v["domain"].map(wtm)
    v = v[v["wtdG"].notna()]
    v["cost"] = v["wtdG"] - v["dG"]
    ss = []
    for dom in sorted(v["domain"].unique()):
        p = os.path.join(afdir, dom.replace(".pdb", "") + ".pdb")
        if not os.path.exists(p):
            continue
        res = D17.parse_backbone(p)
        s, keys = D17.assign(res)
        if not s:
            continue
        el = D17.elements(s, keys)
        for i in keys:
            ss.append({"domain": dom, "pos": i, "aa_model": res[i]["aa"], "ss": el[i][0]})
    d = v.groupby(["domain", "pos", "aa"])["cost"].mean().reset_index().merge(
        pd.DataFrame(ss), on=["domain", "pos"])
    d = d[(d["aa"] == d["aa_model"]) & d["aa"].isin(PACE)].copy()
    return d


def per_residue_relative(d, cls, mode="centre"):
    """Mean cost of deleting each residue from class `cls`, relative to its domain's own mean."""
    cells = defaultdict(list)
    for dom, g in d.groupby("domain"):
        base = float(g["cost"].mean())
        scale = float(g["cost"].std()) if mode == "z" else 1.0
        if mode == "z" and (not np.isfinite(scale) or scale < 1e-6):
            continue
        sub = g[g["ss"] == cls]
        for aa, h in sub.groupby("aa"):
            cells[str(aa)].append((float(h["cost"].mean()) - base) / scale)
    return {a: (float(np.mean(v)), len(v)) for a, v in cells.items() if len(v) >= 5}


def rho(pts):
    aas = [a for a in pts if a in PACE]
    if len(aas) < 8:
        return np.nan, len(aas)
    return float(stats.spearmanr([PACE[a] for a in aas],
                                 [pts[a][0] for a in aas])[0]), len(aas)


def boot_rho(d, cls, mode="centre", nboot=800):
    doms = d["domain"].unique()
    idx = {x: d.index[d["domain"] == x].to_numpy() for x in doms}
    out = []
    for _ in range(nboot):
        pick = RNG.choice(doms, len(doms), replace=True)
        sub = d.loc[np.concatenate([idx[x] for x in pick])]
        r, _ = rho(per_residue_relative(sub, cls, mode))
        if np.isfinite(r):
            out.append(r)
    if not out:
        return [None, None]
    return [round(float(q), 4) for q in np.quantile(out, [.025, .975])]


def main(lib3, afdir, out_path):
    d = load(lib3, afdir)
    rep = {"n": int(len(d)), "domains": int(d["domain"].nunique()),
           "residues": int(d["aa"].nunique()),
           "note": "cost relative to each domain's own mean deletion cost; proline excluded"}
    for mode in ("centre", "z"):
        block = {}
        for cls, nm in (("H", "helix"), ("L", "loop"), ("E", "strand")):
            pts = per_residue_relative(d, cls, mode)
            r, n = rho(pts)
            block[nm] = {"spearman_vs_helix_propensity": round(float(r), 4), "n_residues": n,
                         "ci95_by_domain": boot_rho(d, cls, mode, nboot=500 if mode == "centre"
                                                    else 300),
                         "values": {k: round(v[0], 4) for k, v in sorted(pts.items())}}
        rep["relative_to_domain_mean" if mode == "centre" else "z_scored"] = block

    # the contrast itself, for reference, built from the same cells
    hp = per_residue_relative(d, "H")
    lp = per_residue_relative(d, "L")
    both = sorted(set(hp) & set(lp))
    contrast = {a: hp[a][0] - lp[a][0] for a in both}
    r = (float(stats.spearmanr([PACE[a] for a in both], [contrast[a] for a in both])[0])
         if len(both) >= 8 else np.nan)
    rep["contrast_from_same_cells"] = {"spearman": round(float(r), 4), "n_residues": len(both),
                                       "values": {k: round(v, 4) for k, v in contrast.items()}}

    json.dump(rep, open(out_path, "w"), indent=2)
    print(json.dumps(rep, indent=2))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])
