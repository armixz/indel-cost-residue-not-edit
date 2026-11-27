"""Is the slope of deletion cost on helix propensity the one thermodynamics predicts?

Helix propensity, as measured by Pace & Scholtz, is the free energy a residue contributes to
helix formation relative to alanine: a residue with a penalty of 0.6 kcal/mol stabilises a helix
0.6 kcal/mol less well than alanine does. Deleting a residue from a helix removes whatever that
residue was contributing. If the rest of the helix is otherwise unperturbed, the extra cost of
deleting residue X from a helix rather than from a loop should therefore fall, one for one, as X's
propensity penalty rises:

    cost_helix(X) - cost_loop(X)  =  c - 1.0 * penalty(X)

So this is not only a rank prediction. It fixes the SLOPE at -1 kcal/mol per kcal/mol. A rank
correlation is consistent with any monotone relationship; a slope near -1 would say the effect is
the residue's own helical contribution and nothing else, while a slope far from -1 would say
something additional is going on.

Fitted here by ordinary least squares over residues, weighted by the number of contributing
proteins, with a bootstrap over proteins for the interval. Proline is excluded (not in the
published scale). Reported both with and without glycine, which sits at the end of the range and
could lever the slope on its own.

Usage: indel_env python 45_slope.py <lib3_csv> <afdir> <out_json>
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
    return d[(d["aa"] == d["aa_model"]) & d["aa"].isin(PACE)].copy()


def cells(d, a="H", b="L"):
    out = {}
    for (dom, aa), g in d.groupby(["domain", "aa"]):
        x, y = g[g["ss"] == a]["cost"], g[g["ss"] == b]["cost"]
        if len(x) >= 1 and len(y) >= 1:
            out.setdefault(str(dom), {})[str(aa)] = float(x.mean() - y.mean())
    return out


def fit(cell, picked, drop=(), min_dom=5, weighted=True):
    acc = defaultdict(list)
    for dm in picked:
        for aa, val in cell[dm].items():
            acc[aa].append(val)
    pts = [(PACE[a], float(np.mean(v)), len(v)) for a, v in acc.items()
           if a in PACE and a not in drop and len(v) >= min_dom]
    if len(pts) < 8:
        return None
    x = np.array([p[0] for p in pts])
    y = np.array([p[1] for p in pts])
    w = np.array([p[2] for p in pts], dtype=float) if weighted else np.ones(len(pts))
    W = np.diag(w / w.sum())
    X = np.column_stack([np.ones(len(x)), x])
    beta = np.linalg.solve(X.T @ W @ X, X.T @ W @ y)
    rho = float(stats.spearmanr(x, y)[0])
    return {"intercept": float(beta[0]), "slope": float(beta[1]), "rho": rho,
            "n_residues": len(pts)}


def boot(cell, doms, drop=(), nboot=2000, weighted=True):
    out = defaultdict(list)
    for _ in range(nboot):
        f = fit(cell, list(RNG.choice(doms, len(doms), replace=True)), drop,
                weighted=weighted)
        if f:
            for k in ("intercept", "slope", "rho"):
                out[k].append(f[k])
    return {k: [round(float(q), 4) for q in np.quantile(v, [.025, .975])]
            for k, v in out.items() if v}


def main(lib3, afdir, out_path):
    d = load(lib3, afdir)
    cell = cells(d)
    doms = sorted(cell)
    rep = {"n_positions": int(len(d)), "n_domains_with_contrast": len(doms),
           "prediction": "slope = -1.0 kcal/mol per kcal/mol if the helix contrast is exactly "
                         "the residue's own helical contribution"}
    for lab, drop, w in (("weighted_all", (), True),
                         ("weighted_no_gly", ("G",), True),
                         ("unweighted_all", (), False)):
        f = fit(cell, doms, drop, weighted=w)
        if not f:
            rep[lab] = {"status": "too few residues"}
            continue
        b = boot(cell, doms, drop, weighted=w)
        rep[lab] = {"slope": round(f["slope"], 4), "slope_ci95": b.get("slope"),
                    "intercept": round(f["intercept"], 4),
                    "intercept_ci95": b.get("intercept"),
                    "spearman": round(f["rho"], 4), "spearman_ci95": b.get("rho"),
                    "n_residues": f["n_residues"],
                    "slope_ci_contains_minus_one": bool(
                        b.get("slope") and b["slope"][0] <= -1.0 <= b["slope"][1])}

    # the same for the strand contrast, where no verified propensity scale exists; reported as a
    # descriptive comparison only
    cellE = cells(d, "E", "L")
    domsE = sorted(cellE)
    fE = fit(cellE, domsE)
    if fE:
        bE = boot(cellE, domsE)
        rep["strand_contrast_vs_helix_propensity"] = {
            "slope": round(fE["slope"], 4), "slope_ci95": bE.get("slope"),
            "spearman": round(fE["rho"], 4), "n_residues": fE["n_residues"],
            "note": "helix propensity should NOT predict the strand contrast; reported as a "
                    "negative control on the scale itself"}

    # per-residue table for the manuscript
    acc = defaultdict(list)
    for dm in doms:
        for aa, val in cell[dm].items():
            acc[aa].append(val)
    rep["per_residue"] = {a: [round(float(np.mean(v)), 4), len(v), PACE[a]]
                          for a, v in sorted(acc.items()) if a in PACE and len(v) >= 5}

    json.dump(rep, open(out_path, "w"), indent=2)
    print(json.dumps(rep, indent=2))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])
