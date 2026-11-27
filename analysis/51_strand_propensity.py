"""Does the strand arm have its own residue ordering, matched to beta propensity?

Pre-registered in PREREG_strand_propensity.md, including the addendum recording that the
context precondition is TRIGGERED: because beta propensity is known to be strongly
context-dependent (Minor & Kim, Nature 1994;371:264), a NULL result here is uninformative
and only a positive result is interpretable. That was written before any number below was
computed.

Scales, verified against photographic reproductions of the printed tables, and aligned so
that LARGER = WORSE former, matching the manuscript's helix penalty scale. The two GB1
scales use opposite sign conventions as printed, which is the main trap here:

  Minor & Kim 1994a   Nature 367:660     GB1 residue 53, central strand, Ala = 0
                      printed larger = BETTER, so penalty = -(printed)
  Smith/Withka/Regan  Biochemistry 33:5510  GB1 residue 53, central strand, Ala = 0
                      printed more negative = BETTER, so penalty = +(printed)
  Kim & Berg 1993     Nature 362:267     zinc finger, Gly = 0
                      printed more negative = BETTER, so penalty = +(printed)
                      range 4-5x compressed; a third check, never the primary

A sanity check fixed in the pre-registration runs first: after alignment glycine must
carry the highest penalty on every scale, or the analysis stops.

Usage: python 51_strand_propensity.py <lib3_csv> <afdir> <results_dir> <out_json>
"""
import importlib.util
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
sp = importlib.util.spec_from_file_location("d43", os.path.join(HERE, "43_discriminate.py"))
D43 = importlib.util.module_from_spec(sp)
sp.loader.exec_module(D43)

RNG = np.random.default_rng(0)
NBOOT = 2000
MIN_DOM = 5

PACE = dict(D43.PACE)   # helix penalty, larger = worse

# printed values, before alignment
_MK = {"T": 1.1, "I": 1.0, "Y": 0.96, "F": 0.86, "V": 0.82, "M": 0.72, "S": 0.70,
       "W": 0.54, "C": 0.52, "L": 0.51, "R": 0.45, "K": 0.27, "Q": 0.23, "E": 0.01,
       "A": 0.00, "H": -0.02, "N": -0.08, "D": -0.94, "G": -1.2}
_SR = {"Y": -1.63, "T": -1.36, "I": -1.25, "F": -1.08, "W": -1.04, "V": -0.94,
       "S": -0.87, "M": -0.90, "C": -0.78, "L": -0.45, "R": -0.40, "N": -0.52,
       "H": -0.37, "Q": -0.38, "K": -0.35, "E": -0.23, "A": 0.0, "D": 0.85, "G": 1.21}
_KB = {"I": -0.56, "F": -0.55, "V": -0.53, "Y": -0.50, "W": -0.48, "T": -0.48,
       "L": -0.48, "C": -0.47, "M": -0.46, "H": -0.46, "R": -0.44, "K": -0.41,
       "E": -0.41, "D": -0.41, "Q": -0.40, "S": -0.39, "N": -0.38, "A": -0.35, "G": 0.0}

BETA = {"minor_kim_central": {a: -v for a, v in _MK.items()},
        "smith_regan_central": dict(_SR),
        "kim_berg_zincfinger": dict(_KB)}


def check_alignment():
    bad = []
    for name, sc in BETA.items():
        if max(sc, key=sc.get) != "G":
            bad.append(f"{name}: highest penalty is {max(sc, key=sc.get)}, not G")
    return bad


def cells(d, a, b, extra=None):
    """{(domain, aa): mean cost in class a minus mean cost in class b}, paired in domain."""
    out = {}
    src = d if extra is None else d[extra]
    for (dom, aa), g in src.groupby(["domain", "aa"]):
        x, y = g[g["cls"] == a]["cost"], g[g["cls"] == b]["cost"]
        if len(x) >= 1 and len(y) >= 1:
            out[(str(dom), str(aa))] = float(x.mean() - y.mean())
    return out


def per_aa(cell, doms):
    acc = {}
    keep = set(doms)
    for (dom, aa), v in cell.items():
        if dom in keep:
            acc.setdefault(aa, []).append(v)
    return {aa: float(np.mean(v)) for aa, v in acc.items() if len(v) >= MIN_DOM}


def rho(contrasts, scale):
    k = sorted(a for a in contrasts if a in scale)
    if len(k) < 8:
        return np.nan
    return float(stats.spearmanr([scale[a] for a in k], [contrasts[a] for a in k])[0])


def ci(vals):
    vals = [v for v in vals if np.isfinite(v)]
    if not vals:
        return [None, None]
    return [round(float(q), 4) for q in np.quantile(vals, [.025, .975])]


def main(lib3, afdir, rdir, out_path):
    bad = check_alignment()
    rep = {"prereg": "PREREG_strand_propensity.md",
           "context_precondition": "TRIGGERED before any computation: a null result on S1 "
                                   "or S2 is uninformative, because beta propensity is "
                                   "known to be strongly context dependent (Minor & Kim, "
                                   "Nature 1994;371:264). Only a positive result is "
                                   "interpretable.",
           "alignment_sanity_check": "PASS" if not bad else bad}
    if bad:
        json.dump(rep, open(out_path, "w"), indent=2)
        print("ALIGNMENT CHECK FAILED:", bad)
        return

    d = D43.load(lib3, afdir, rdir)
    d["cls"] = d["ss"]

    # edge vs central strands, from bridge partner counts
    br = pd.read_csv(os.path.join(rdir, "tsub_bridges.csv"))
    br["domain"] = br["domain"].astype(str) + ".pdb"
    br = br.rename(columns={"struct_pos": "pos"})[["domain", "pos", "n_partners"]]
    d = d.merge(br, on=["domain", "pos"], how="left")
    d["strand_kind"] = np.where(d["cls"] != "E", "",
                                np.where(d["n_partners"] >= 2, "central",
                                         np.where(d["n_partners"] == 1, "edge", "unknown")))
    rep["strand_positions"] = {
        k: int(v) for k, v in d[d["cls"] == "E"]["strand_kind"].value_counts().items()}

    doms = sorted(d["domain"].unique())
    cE = cells(d, "E", "L")
    cH = cells(d, "H", "L")
    dEc = d[(d["cls"] == "L") | (d["strand_kind"] == "central")]
    dEe = d[(d["cls"] == "L") | (d["strand_kind"] == "edge")]
    cEc = cells(dEc, "E", "L")
    cEe = cells(dEe, "E", "L")

    scales = {"pace_helix": PACE}
    scales.update(BETA)

    point = {}
    for cname, cell in (("strand", cE), ("helix", cH),
                        ("strand_central", cEc), ("strand_edge", cEe)):
        ctr = per_aa(cell, doms)
        point[cname] = {"n_residues": len(ctr),
                        **{s: round(rho(ctr, sc), 4) for s, sc in scales.items()}}

    boots = {c: {s: [] for s in scales} for c in point}
    diffs = {"S2_strand_beta_minus_helixscale": [], "S3_helix_helixscale_minus_beta": [],
             "S6_central_minus_edge_beta": []}
    PRIMARY = "minor_kim_central"
    for _ in range(NBOOT):
        pick = list(RNG.choice(doms, len(doms), replace=True))
        ctrs = {c: per_aa(cell, pick) for c, cell in
                (("strand", cE), ("helix", cH), ("strand_central", cEc), ("strand_edge", cEe))}
        rr = {}
        for c, ctr in ctrs.items():
            for s, sc in scales.items():
                v = rho(ctr, sc)
                boots[c][s].append(v)
                rr[(c, s)] = v
        # magnitudes: a stronger relationship is a more negative rho
        diffs["S2_strand_beta_minus_helixscale"].append(
            rr[("strand", PRIMARY)] - rr[("strand", "pace_helix")])
        diffs["S3_helix_helixscale_minus_beta"].append(
            rr[("helix", "pace_helix")] - rr[("helix", PRIMARY)])
        diffs["S6_central_minus_edge_beta"].append(
            rr[("strand_central", PRIMARY)] - rr[("strand_edge", PRIMARY)])

    rep["S1_S4_correlations"] = {
        c: {"n_residues": point[c]["n_residues"],
            **{s: {"spearman": point[c][s], "ci95": ci(boots[c][s])} for s in scales}}
        for c in point}

    rep["S2_S3_S6_differences"] = {
        k: {"mean": round(float(np.nanmean(v)), 4), "ci95": ci(v),
            "excludes_zero": bool(ci(v)[0] is not None
                                  and (ci(v)[0] > 0 or ci(v)[1] < 0))}
        for k, v in diffs.items()}
    rep["difference_sign_note"] = (
        "A more negative rho means a stronger relationship. S2 is negative when beta "
        "propensity beats helix propensity on the strand arm; S3 is negative when helix "
        "propensity beats beta propensity on the helix arm; S6 is negative when the beta "
        "relationship is stronger at central strands than at edge strands, which is the "
        "pre-registered prediction.")

    # S5: how collinear are the scales themselves?
    rep["S5_scale_collinearity"] = {}
    for s, sc in BETA.items():
        k = sorted(set(sc) & set(PACE))
        rep["S5_scale_collinearity"][f"pace_vs_{s}"] = {
            "n_residues": len(k),
            "spearman": round(float(stats.spearmanr([PACE[a] for a in k],
                                                    [sc[a] for a in k])[0]), 4)}
    k = sorted(set(BETA["minor_kim_central"]) & set(BETA["smith_regan_central"]))
    rep["S5_scale_collinearity"]["minor_kim_vs_smith_regan"] = {
        "n_residues": len(k),
        "spearman": round(float(stats.spearmanr(
            [BETA["minor_kim_central"][a] for a in k],
            [BETA["smith_regan_central"][a] for a in k])[0]), 4)}

    json.dump(rep, open(out_path, "w"), indent=2)

    print(f"alignment check: {rep['alignment_sanity_check']}")
    print(f"strand positions: {rep['strand_positions']}\n")
    print(f"{'contrast':16s} {'n':>3s} " + " ".join(f"{s[:13]:>14s}" for s in scales))
    for c in point:
        print(f"{c:16s} {point[c]['n_residues']:3d} " +
              " ".join(f"{point[c][s]:14.3f}" for s in scales))
    print("\ndifferences (negative = the first term is the stronger relationship):")
    for k_, v in rep["S2_S3_S6_differences"].items():
        print(f"  {k_:36s} {v['mean']:7.3f}  {v['ci95']}  excludes zero: {v['excludes_zero']}")
    print("\nscale collinearity (Spearman between the scales themselves):")
    for k_, v in rep["S5_scale_collinearity"].items():
        print(f"  {k_:34s} {v['spearman']:6.3f}  (n={v['n_residues']})")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
