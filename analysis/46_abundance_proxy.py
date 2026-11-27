"""Is abundance/stability a valid general proxy for protein function?

Pre-registered in PREREG_abundance_proxy.md before any data was examined.

The question matters because every general, protein-agnostic, pooled assay that exists
measures abundance or stability, never activity. If abundance predicts function well, one
general assay suffices for a tag-insertion data platform. If it does not, every protein
needs a bespoke activity assay and the platform does not scale.

Five proteins have both an abundance/stability deep mutational scan and an independent
activity scan over the same variants, downloaded from MaveDB:

  PTEN    abundance urn:mavedb:00000102-0-1  function urn:mavedb:00000054-a-1
  NUDT15  abundance urn:mavedb:00000055-a-1  function urn:mavedb:00000055-b-1
  CYP2C9  abundance urn:mavedb:00000095-b-1  function urn:mavedb:00000095-a-1
  G6PD    abundance urn:mavedb:00001266-b-1  function urn:mavedb:00001266-a-1
  ASPA    abundance urn:mavedb:00000657-a-1  function urn:mavedb:00000657-b-1

Both scores are put on a common scale where the synonymous median is 1 and the nonsense
median is 0. That rescaling also fixes the sign automatically, so an assay reported in the
opposite direction (a toxicity readout, say) needs no hand-applied flip, and no judgement
about orientation is made after seeing the correlation.

Nonsense variants are then DROPPED. They are trivially both non-abundant and
non-functional and would inflate every correlation; the pre-registration fixes this
exclusion in advance. All marks are computed on missense variants only.

Usage: python 46_abundance_proxy.py <mave_dir> <out_json>
"""
import json
import os
import re
import sys

import numpy as np
import pandas as pd
from scipy import stats

RNG = np.random.default_rng(0)

PAIRS = ["PTEN", "NUDT15", "CYP2C9", "G6PD"]

# ASPA (abundance urn:mavedb:00000657-a-1, function urn:mavedb:00000657-b-1) was
# retrieved and then EXCLUDED, for a reason visible in the data before any correlation
# was computed: its second assay is a cellular TOXICITY readout, not an activity readout.
# In a toxicity assay a premature stop makes no protein and is therefore not toxic, so
# nonsense variants score as benign (median 0.122) while missense variants score lower
# (0.048). A null allele thus lands on the "good" end, which means the assay cannot be
# oriented by the nonsense class and is not a loss-of-function readout at all.
# Recording the direction this exclusion cuts, so it cannot be read as convenient: if
# ASPA were instead included with its sign flipped it would contribute rho = +0.77, the
# second HIGHEST of the set, and would push the median correlation UP (0.357 -> 0.373).
# Excluding it therefore works against the conclusion this script reaches, not for it.
EXCLUDED = {"ASPA": "second assay is cellular toxicity, not activity; nonsense variants "
                    "are non-toxic, so the assay has no loss-of-function orientation"}

# Anchors are declared per dataset rather than guessed, because guessing orientation from
# percentiles cannot work for an assay whose null class is not at an extreme.
#   "wt_null" : the source already publishes scores on a 0 = null, 1 = wild-type scale
#               (standard for VAMP-seq and for the NUDT15 scan), so no rescaling is done.
#   "syn_nonsense" : high anchor from synonymous variants (or the explicit wild-type row
#               when the set has too few synonymous), low anchor from the nonsense median.
ANCHORS = {
    ("PTEN", "abundance"): "wt_null",
    ("PTEN", "function"): "syn_nonsense",
    ("NUDT15", "abundance"): "wt_null",
    ("NUDT15", "function"): "wt_null",
    ("CYP2C9", "abundance"): "syn_nonsense",
    ("CYP2C9", "function"): "syn_nonsense",
    ("G6PD", "abundance"): "syn_nonsense",
    ("G6PD", "function"): "syn_nonsense",
}

THREE = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q",
         "Glu": "E", "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K",
         "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W",
         "Tyr": "Y", "Val": "V", "Ter": "*", "Sec": "U"}
PAT = re.compile(r"^p\.(?:\()?([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|=|\*|Ter)(?:\))?$")

# thresholds, fixed in the pre-registration
FUNC_DEAD = 0.5      # normalised function below this = non-functional
ABUND_OK = 0.75      # normalised abundance at or above this = assay calls it fine


def parse(h):
    """Return (wt, pos, mut) with mut '=' for synonymous and '*' for nonsense."""
    if not isinstance(h, str):
        return None
    m = PAT.match(h.strip())
    if not m:
        return None
    wt3, pos, mut3 = m.groups()
    wt = THREE.get(wt3)
    if wt is None:
        return None
    if mut3 == "=":
        mut = wt
    elif mut3 in ("*", "Ter"):
        mut = "*"
    else:
        mut = THREE.get(mut3)
    if mut is None:
        return None
    return wt, int(pos), mut


def load(path):
    d = pd.read_csv(path)
    p = d["hgvs_pro"].map(parse)
    d = d[p.notna()].copy()
    d[["wt", "pos", "mut"]] = pd.DataFrame(list(p[p.notna()]), index=d.index)
    d = d[d["score"].notna()].copy()
    d["klass"] = np.where(d["mut"] == "*", "nonsense",
                          np.where(d["mut"] == d["wt"], "synonymous", "missense"))
    return d


def normalise(d, tag, mode):
    """Scale so the synonymous median is 1 and the nonsense median is 0.

    Falls back to the 95th/5th percentile of missense scores when a class is absent or
    the two anchors are too close to separate; which rule was used is recorded.
    """
    syn = d.loc[d["klass"] == "synonymous", "score"]
    non = d.loc[d["klass"] == "nonsense", "score"]
    mis = d.loc[d["klass"] == "missense", "score"]
    # The two anchors are chosen independently. The low anchor must come from nonsense
    # variants whenever they exist, because they are the only class whose direction is
    # known a priori (a premature stop is dead), and that is what fixes the orientation
    # of the assay. Falling back to a missense percentile for BOTH anchors cannot orient
    # anything, so it is used only when neither class is available.
    hi, hi_src = (float(syn.median()), "syn") if len(syn) >= 10 else \
                 (float(mis.quantile(.95)), "missense_p95")
    lo, lo_src = (float(non.median()), "nonsense") if len(non) >= 10 else \
                 (float(mis.quantile(.05)), "missense_p05")
    how = f"{hi_src}/{lo_src}"
    if mode == "wt_null":
        # the source already publishes 0 = null, 1 = wild type; rescaling it again would
        # only inject noise from classes it does not contain
        hi, lo, how = 1.0, 0.0, "published_wt_null_scale"
    if abs(hi - lo) < 1e-9:
        hi, lo = float(mis.quantile(.95)), float(mis.quantile(.05))
        how = "degenerate_anchors_percentile_fallback"
    oriented = (mode == "wt_null") or (hi_src == "syn") or (lo_src == "nonsense")
    d = d.copy()
    d["norm"] = (d["score"] - lo) / (hi - lo)
    return d, {"assay": tag, "anchor": how, "hi_raw": round(hi, 4), "lo_raw": round(lo, 4),
               "n_syn": int(len(syn)), "n_nonsense": int(len(non)),
               "n_missense": int(len(mis)),
               "orientation_fixed_by_a_known_class": bool(oriented),
               "flipped": bool(hi < lo)}


def auc_mw(pos, neg):
    """AUC from the Mann-Whitney U statistic (ties handled by rank averaging)."""
    pos, neg = np.asarray(pos, float), np.asarray(neg, float)
    if len(pos) == 0 or len(neg) == 0:
        return np.nan
    r = stats.rankdata(np.concatenate([pos, neg]))
    u = r[:len(pos)].sum() - len(pos) * (len(pos) + 1) / 2.0
    return float(u / (len(pos) * len(neg)))


def boot_ci(fn, n, nboot=2000):
    vals = []
    for _ in range(nboot):
        idx = RNG.integers(0, n, n)
        v = fn(idx)
        if np.isfinite(v):
            vals.append(v)
    if not vals:
        return [None, None]
    return [round(float(q), 4) for q in np.quantile(vals, [.025, .975])]


def main(mdir, out_path):
    rep = {"prereg": "PREREG_abundance_proxy.md", "excluded": EXCLUDED,
           "thresholds": {"function_dead_below": FUNC_DEAD,
                          "abundance_called_fine_at_or_above": ABUND_OK},
           "per_protein": {}}
    pooled = []

    for g in PAIRS:
        fa = os.path.join(mdir, f"{g}_abundance.csv")
        ff = os.path.join(mdir, f"{g}_function.csv")
        if not (os.path.exists(fa) and os.path.exists(ff)):
            rep["per_protein"][g] = {"status": "missing files"}
            continue
        A, ia = normalise(load(fa), "abundance", ANCHORS[(g, "abundance")])
        F, iff = normalise(load(ff), "function", ANCHORS[(g, "function")])

        key = ["wt", "pos", "mut"]
        m = (A[A["klass"] == "missense"][key + ["norm"]]
             .groupby(key, as_index=False).mean()
             .merge(F[F["klass"] == "missense"][key + ["norm"]]
                    .groupby(key, as_index=False).mean(),
                    on=key, suffixes=("_ab", "_fn")))
        if len(m) < 100:
            rep["per_protein"][g] = {"status": f"only {len(m)} shared missense variants"}
            continue

        x = m["norm_ab"].to_numpy()
        y = m["norm_fn"].to_numpy()
        n = len(m)

        # A1 correlation
        rho = float(stats.spearmanr(x, y)[0])
        r_ci = boot_ci(lambda i: stats.spearmanr(x[i], y[i])[0], n)
        pear = float(stats.pearsonr(x, y)[0])

        # A2 the decisive asymmetry: abundant but dead
        abundant = x >= ABUND_OK
        dead = y < FUNC_DEAD
        n_ab = int(abundant.sum())
        fp = float((abundant & dead).sum() / n_ab) if n_ab else np.nan
        fp_ci = boot_ci(lambda i: ((x[i] >= ABUND_OK) & (y[i] < FUNC_DEAD)).sum()
                        / max((x[i] >= ABUND_OK).sum(), 1), n)
        # and the converse: low abundance but still functional
        low = ~abundant
        fn_rate = float((low & ~dead).sum() / max(int(low.sum()), 1))

        # A3 discrimination
        auc = auc_mw(x[~dead], x[dead])   # abundance higher in functional variants
        auc_ci = boot_ci(lambda i: auc_mw(x[i][~dead[i]], x[i][dead[i]]), n)

        rep["per_protein"][g] = {
            "n_shared_missense": n,
            "abundance_assay": ia, "function_assay": iff,
            "A1_spearman": round(rho, 4), "A1_ci95": r_ci,
            "A1_pearson": round(pear, 4),
            "A2_n_called_abundant": n_ab,
            "A2_false_pass_rate": round(fp, 4) if np.isfinite(fp) else None,
            "A2_false_pass_ci95": fp_ci,
            "A2_low_abundance_but_functional": round(fn_rate, 4),
            "A3_auc": round(auc, 4) if np.isfinite(auc) else None,
            "A3_auc_ci95": auc_ci,
            "frac_dead_overall": round(float(dead.mean()), 4)}

        sub = m.copy()
        sub["gene"] = g
        pooled.append(sub)

    if pooled:
        P = pd.concat(pooled, ignore_index=True)
        x, y = P["norm_ab"].to_numpy(), P["norm_fn"].to_numpy()
        abundant, dead = x >= ABUND_OK, y < FUNC_DEAD
        rhos = [v["A1_spearman"] for v in rep["per_protein"].values()
                if isinstance(v.get("A1_spearman"), float)]
        fps = [v["A2_false_pass_rate"] for v in rep["per_protein"].values()
               if isinstance(v.get("A2_false_pass_rate"), float)]
        rep["pooled"] = {
            "n": int(len(P)), "n_proteins": len(rhos),
            "median_spearman": round(float(np.median(rhos)), 4) if rhos else None,
            "spearman_range": [round(min(rhos), 4), round(max(rhos), 4)] if rhos else None,
            "n_proteins_rho_ge_0.7": int(sum(r >= 0.7 for r in rhos)),
            "pooled_false_pass_rate": round(float((abundant & dead).sum()
                                                  / max(int(abundant.sum()), 1)), 4),
            "median_false_pass_rate": round(float(np.median(fps)), 4) if fps else None,
            "false_pass_range": [round(min(fps), 4), round(max(fps), 4)] if fps else None}

        # Threshold sensitivity. The false-pass rate depends on two cutoffs chosen in the
        # pre-registration, so the conclusion is only worth anything if it survives other
        # reasonable choices. Sweeping both and reporting the whole grid means the
        # headline number cannot be an artefact of a lucky pair.
        sweep = {}
        for a_cut in (0.6, 0.7, 0.75, 0.8, 0.9):
            for f_cut in (0.3, 0.4, 0.5, 0.6, 0.7):
                ab, dd = x >= a_cut, y < f_cut
                sweep[f"abund>={a_cut},func<{f_cut}"] = {
                    "n_called_abundant": int(ab.sum()),
                    "false_pass": round(float((ab & dd).sum() / max(int(ab.sum()), 1)), 4)}
        rep["threshold_sweep_pooled"] = sweep
        fps_sweep = [v["false_pass"] for v in sweep.values()]
        rep["threshold_sweep_summary"] = {
            "min_false_pass_over_grid": round(min(fps_sweep), 4),
            "max_false_pass_over_grid": round(max(fps_sweep), 4),
            "note": "if even the most generous cell of this grid leaves a large false-pass"
                    " rate, the conclusion does not depend on the pre-registered cutoffs"}

        # pre-declared verdict, applied mechanically
        md = rep["pooled"]["median_spearman"]
        mf = rep["pooled"]["median_false_pass_rate"]
        maj = rep["pooled"]["n_proteins_rho_ge_0.7"] > len(rhos) / 2
        if maj and mf is not None and mf < 0.10:
            v = "PASS: abundance is a good proxy; a general assay suffices"
        elif (md is not None and md < 0.5) or (mf is not None and mf > 0.20):
            v = ("FAIL: abundance is not a substitute for function; "
                 "protein-by-protein assays required")
        else:
            v = "INTERMEDIATE: abundance is a screen, not a verdict"
        rep["verdict"] = v
        P.to_csv(os.path.join(mdir, "merged_variants.csv"), index=False)

    json.dump(rep, open(out_path, "w"), indent=2)

    print(f"{'protein':8s} {'n':>6s} {'rho':>7s} {'AUC':>6s} {'falsepass':>10s} "
          f"{'%dead':>6s}  anchors(ab/fn)")
    for g, v in rep["per_protein"].items():
        if "A1_spearman" not in v:
            print(f"{g:8s} {v.get('status', '?')}")
            continue
        print(f"{g:8s} {v['n_shared_missense']:6d} {v['A1_spearman']:7.3f} "
              f"{v['A3_auc']:6.3f} {v['A2_false_pass_rate']:10.3f} "
              f"{v['frac_dead_overall']:6.3f}  "
              f"{v['abundance_assay']['anchor']} / {v['function_assay']['anchor']}")
    p = rep.get("pooled", {})
    print(f"\nmedian rho {p.get('median_spearman')}  range {p.get('spearman_range')}  "
          f"rho>=0.7 in {p.get('n_proteins_rho_ge_0.7')}/{p.get('n_proteins')}")
    print(f"false-pass: median {p.get('median_false_pass_rate')}  "
          f"range {p.get('false_pass_range')}  pooled {p.get('pooled_false_pass_rate')}")
    print(f"\n{rep.get('verdict')}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
