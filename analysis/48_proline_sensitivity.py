"""Does the propensity result depend on proline, whose scale value we cannot verify?

The manuscript uses the 19 helix propensities listed in the abstract of Pace & Scholtz
(1998), all of which were checked against that abstract. Proline is not among them. The
abstract says "excluding proline, glycine has the lowest helix propensity", which is
ambiguous: it may mean proline is absent from the scale, or that it is present with a
value worse than glycine's. A value of 3.16 kcal/mol circulates in secondary sources but
could not be verified against the paper, whose full text is behind an anti-automation
interstitial on PMC.

Rather than depend on a number we cannot check, this script asks whether it matters.
Proline IS present in the deletion data -- natural domains contain prolines and the
library deletes them -- so its contrast can be computed from our own measurements. The
test is then:

  1. the headline Spearman over the 19 verified residues (the manuscript's number)
  2. the same with proline added at the unverified 3.16
  3. the same with proline added at a range of hypothetical values, to see whether ANY
     plausible proline value would change the conclusion
  4. where proline's measured contrast actually falls, which is itself informative: if
     proline behaves like a poor helix former, that is a consistency check on the whole
     propensity story, obtained without needing the published number at all

Usage: python 48_proline_sensitivity.py <lib3_csv> <afdir> <results_dir> <out_json>
"""
import importlib.util
import json
import os
import sys

import numpy as np
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
sp = importlib.util.spec_from_file_location("d43", os.path.join(HERE, "43_discriminate.py"))
D43 = importlib.util.module_from_spec(sp)
sp.loader.exec_module(D43)

PACE19 = dict(D43.PACE)          # the verified 19, proline absent
PRO_CIRCULATED = 3.16


def main(lib3, afdir, rdir, out_path):
    rep = {}

    # The manuscript's own frame: D43.load drops residues absent from the Pace scale,
    # which removes proline.
    d = D43.load(lib3, afdir, rdir)
    pts = D43.contrasts(d)
    aas = sorted(pts)
    y = np.array([pts[a][0] for a in aas])
    x19 = np.array([PACE19[a] for a in aas])
    rho19 = float(stats.spearmanr(x19, y)[0])
    rep["verified_19"] = {"n_residues": len(aas), "spearman": round(rho19, 4),
                          "residues": aas}

    # The same frame with proline kept, so its own contrast can be measured from our
    # data without needing any published value for it.
    full = _load_with_proline(lib3, afdir, rdir)
    pts_full = D43.contrasts(full)
    rep["proline_measured"] = (
        {"contrast_kcal": round(pts_full["P"][0], 4), "n_domains": pts_full["P"][1]}
        if "P" in pts_full else {"status": "too few domains with both P-in-helix and P-in-loop"})
    gly = pts_full.get("G")
    rep["glycine_measured"] = ({"contrast_kcal": round(gly[0], 4), "n_domains": gly[1]}
                               if gly else None)
    rep["rank_of_proline_among_measured"] = None
    if "P" in pts_full:
        vals = sorted(((v[0], k) for k, v in pts_full.items()))
        rep["rank_of_proline_among_measured"] = {
            "rank_from_most_negative": [k for _, k in vals].index("P") + 1,
            "n_residues_ranked": len(vals),
            "note": "a poor helix former should sit toward the MOST NEGATIVE end, like "
                    "glycine, because deleting it from a helix costs least"}

    # Sensitivity: add proline at a range of hypothetical scale values.
    if "P" in pts_full:
        aas_p = sorted(pts_full)
        sweep = {}
        for pro_val in (1.00, 1.50, 2.00, 2.50, 3.16, 4.00, 5.00):
            scale = dict(PACE19)
            scale["P"] = pro_val
            keep = [a for a in aas_p if a in scale]
            xx = np.array([scale[a] for a in keep])
            yy = np.array([pts_full[a][0] for a in keep])
            sweep[f"proline={pro_val:.2f}"] = {
                "n_residues": len(keep),
                "spearman": round(float(stats.spearmanr(xx, yy)[0]), 4)}
        rep["sensitivity_sweep"] = sweep
        rhos = [v["spearman"] for v in sweep.values()]
        rep["sensitivity_summary"] = {
            "spearman_without_proline": round(rho19, 4),
            "range_with_proline": [round(min(rhos), 4), round(max(rhos), 4)],
            "max_absolute_change": round(max(abs(r - rho19) for r in rhos), 4)}

    json.dump(rep, open(out_path, "w"), indent=2)
    print(json.dumps(rep, indent=2))


def _load_with_proline(lib3, afdir, rdir):
    """Re-run the manuscript's loader with proline present.

    D43.load ends with d[d['pace'].notna()], so a residue missing from the scale is
    dropped. Temporarily giving proline an entry keeps it through that filter; the entry
    is removed again afterwards so nothing else in this process sees a modified scale.
    The value used here is irrelevant to the measured contrast -- it only has to exist.
    """
    D43.PACE["P"] = PRO_CIRCULATED
    try:
        return D43.load(lib3, afdir, rdir)
    finally:
        D43.PACE.pop("P", None)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
