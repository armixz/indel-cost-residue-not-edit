"""Propensity or size? A discriminating test against the periodicity account.

Goldman, Chati & Ntranos (bioRxiv 2024.06.27.601077) explain the special sensitivity of helices
to indels by helical periodicity: a size-altering edit breaks the register, a size-conservative
one does not, and they report that "non-size-conservative substitutions tend to be more damaging
within helices". On that account the helix penalty is a property of the EDIT (its effect on
register and size), and any apparent residue dependence should track residue SIZE.

The account proposed here is that the helix penalty tracks how well the residue fits a helix,
i.e. its helix propensity. The two make different predictions for which scale explains the
per-residue helix-minus-loop contrast in deletion cost.

Tests:
  1. Spearman of the contrast against helix propensity, residue volume, and hydrophobicity
  2. partial correlations: propensity controlling for volume, and volume controlling for propensity
  3. the same within tertiles of relative solvent accessibility (buried faces versus exposed)
  4. a direct regression of deletion cost on propensity x structure and volume x structure
     together, so the two interactions compete in one model
  5. proline is excluded throughout: the Pace & Scholtz abstract excludes it and the circulating
     value of 3.16 could not be verified against the paper

Usage: indel_env python 43_discriminate.py <lib3_csv> <afdir> <results_dir> <out_json>
"""
import csv
import importlib.util
import json
import os
import re
import sys

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
VOL = {"G": 60.1, "A": 88.6, "S": 89.0, "C": 108.5, "D": 111.1, "N": 114.1,
       "T": 116.1, "E": 138.4, "V": 140.0, "Q": 143.8, "H": 153.2, "M": 162.9, "I": 166.7,
       "L": 166.7, "K": 168.6, "R": 173.4, "F": 189.9, "Y": 193.6, "W": 227.8}
# Kyte-Doolittle hydropathy
KD = {"I": 4.5, "V": 4.2, "L": 3.8, "F": 2.8, "C": 2.5, "M": 1.9, "A": 1.8, "G": -0.4,
      "T": -0.7, "S": -0.8, "W": -0.9, "Y": -1.3, "H": -3.2, "Q": -3.5, "N": -3.5,
      "E": -3.5, "D": -3.5, "K": -3.9, "R": -4.5}


def load(lib3, afdir, rdir):
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
    S = pd.DataFrame(ss)
    d = v.groupby(["domain", "pos", "aa"])["cost"].mean().reset_index().merge(
        S, on=["domain", "pos"])
    d = d[d["aa"] == d["aa_model"]]
    feat = pd.read_csv(f"{rdir}/tsub_features.csv")[["domain", "pos", "rsa", "plddt"]]
    d = d.merge(feat, on=["domain", "pos"], how="left")
    for nm, sc in (("pace", PACE), ("vol", VOL), ("kd", KD)):
        d[nm] = d["aa"].map(sc)
    return d[d["pace"].notna()].copy()


def contrasts(df, a="H", b="L", nmin=1, min_dom=5):
    out = {}
    for aa, g in df.groupby("aa"):
        v = [float(x[x["ss"] == a]["cost"].mean() - x[x["ss"] == b]["cost"].mean())
             for _, x in g.groupby("domain")
             if (x["ss"] == a).sum() >= nmin and (x["ss"] == b).sum() >= nmin]
        if len(v) >= min_dom:
            out[str(aa)] = (float(np.mean(v)), len(v))
    return out


def partial_spearman(x, y, z):
    """Spearman partial correlation of x and y given z, via ranks and residuals."""
    rx, ry, rz = (stats.rankdata(v) for v in (x, y, z))
    bx = np.polyfit(rz, rx, 1)
    by = np.polyfit(rz, ry, 1)
    ex = rx - (bx[0] * rz + bx[1])
    ey = ry - (by[0] * rz + by[1])
    return float(stats.pearsonr(ex, ey)[0])


def main(lib3, afdir, rdir, out_path):
    d = load(lib3, afdir, rdir)
    rep = {"n": int(len(d)), "domains": int(d["domain"].nunique()),
           "residues": int(d["aa"].nunique()),
           "note": "proline excluded throughout: not in the verified Pace & Scholtz list"}

    pts = contrasts(d)
    aas = sorted(pts)
    y = np.array([pts[a][0] for a in aas])
    xs = {"pace": np.array([PACE[a] for a in aas]),
          "vol": np.array([VOL[a] for a in aas]),
          "kd": np.array([KD[a] for a in aas])}
    rep["n_residues_in_contrast"] = len(aas)
    rep["scale_intercorrelation"] = {
        "pace_vs_vol": round(float(stats.spearmanr(xs["pace"], xs["vol"])[0]), 3),
        "pace_vs_kd": round(float(stats.spearmanr(xs["pace"], xs["kd"])[0]), 3),
        "vol_vs_kd": round(float(stats.spearmanr(xs["vol"], xs["kd"])[0]), 3)}

    # 1. marginal correlations, with a bootstrap over domains.
    # One contrast per (domain, residue) is computed once; resampling domains then only
    # re-averages those numbers, which is identical to recomputing and far cheaper.
    cell = {}
    for (dom, aa), g in d.groupby(["domain", "aa"]):
        h, lp = g[g["ss"] == "H"]["cost"], g[g["ss"] == "L"]["cost"]
        if len(h) >= 1 and len(lp) >= 1:
            cell.setdefault(str(dom), {})[str(aa)] = float(h.mean() - lp.mean())
    doms = sorted(cell)

    def rho_cells(picked, scale, min_dom=5):
        acc = {}
        for dm in picked:
            for aa, val in cell[dm].items():
                acc.setdefault(aa, []).append(val)
        pairs = [(scale[a], float(np.mean(v))) for a, v in acc.items()
                 if a in scale and len(v) >= min_dom]
        if len(pairs) < 8:
            return np.nan
        xx, yy = zip(*pairs)
        return float(stats.spearmanr(xx, yy)[0])

    def boot_rho(scale, nboot=2000):
        vals = []
        for _ in range(nboot):
            r = rho_cells(list(RNG.choice(doms, len(doms), replace=True)), scale)
            if np.isfinite(r):
                vals.append(r)
        if not vals:
            return [None, None]
        return [round(float(q), 4) for q in np.quantile(vals, [.025, .975])]

    rep["marginal"] = {}
    for nm, sc in (("helix_propensity", PACE), ("volume", VOL), ("hydropathy", KD)):
        arr = np.array([sc[a] for a in aas])
        rep["marginal"][nm] = {"spearman": round(float(stats.spearmanr(arr, y)[0]), 4),
                               "ci95_by_domain": boot_rho(sc)}

    # 2. partial correlations
    rep["partial"] = {
        "propensity_given_volume": round(partial_spearman(xs["pace"], y, xs["vol"]), 4),
        "volume_given_propensity": round(partial_spearman(xs["vol"], y, xs["pace"]), 4),
        "propensity_given_hydropathy": round(partial_spearman(xs["pace"], y, xs["kd"]), 4),
        "hydropathy_given_propensity": round(partial_spearman(xs["kd"], y, xs["pace"]), 4)}

    # 3. within burial tertiles
    dd = d[d["rsa"].notna()].copy()
    dd["t"] = pd.qcut(dd["rsa"], 3, duplicates="drop")
    strata = {}
    for t, g in dd.groupby("t", observed=True):
        p = contrasts(g, min_dom=4)
        a2 = sorted(p)
        if len(a2) >= 8:
            strata[str(t)] = {
                "n_residues": len(a2),
                "propensity": round(float(stats.spearmanr([PACE[a] for a in a2],
                                                          [p[a][0] for a in a2])[0]), 4),
                "volume": round(float(stats.spearmanr([VOL[a] for a in a2],
                                                      [p[a][0] for a in a2])[0]), 4)}
    rep["by_burial"] = strata

    # 4. both interactions in one regression
    try:
        import statsmodels.formula.api as smf
        r = d.dropna(subset=["rsa"]).copy()
        r["dom"] = r["domain"].astype("category")
        r["isH"] = (r["ss"] == "H").astype(float)
        r["isE"] = (r["ss"] == "E").astype(float)
        for c in ("pace", "vol", "kd", "rsa"):
            r[c] = (r[c] - r[c].mean()) / r[c].std()
        cl = {"groups": r["dom"].cat.codes}
        m = smf.ols("cost ~ C(dom) + isH + isE + pace + vol + kd + rsa"
                    " + pace:isH + vol:isH + kd:isH + rsa:isH"
                    " + pace:isE + vol:isE + kd:isE + rsa:isE", data=r).fit(
            cov_type="cluster", cov_kwds=cl)
        keep = [k for k in m.params.index if ":" in k]
        rep["joint_regression"] = {
            "n": int(len(r)),
            "coef": {k: round(float(m.params[k]), 4) for k in keep},
            "ci95": {k: [round(float(a), 4), round(float(b), 4)]
                     for k, (a, b) in m.conf_int().loc[keep].iterrows()},
            "note": "all scales standardised; cost in kcal/mol; errors clustered by domain. "
                    "A positive pace:isH means a worse helix former is costlier to delete from "
                    "a helix; the sign predicted by the propensity account is NEGATIVE."}
    except ImportError:
        rep["joint_regression"] = {"status": "statsmodels unavailable"}

    rep["contrasts"] = {k: [round(v[0], 4), v[1]] for k, v in sorted(pts.items())}
    json.dump(rep, open(out_path, "w"), indent=2)
    print(json.dumps(rep, indent=2))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
