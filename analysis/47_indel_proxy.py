"""The same question as 46, asked directly of indels instead of substitutions.

Script 46 tests whether abundance predicts function using single amino-acid substitutions,
because that is where paired abundance/activity scans are plentiful. Substitutions are a
weaker perturbation than inserting a peptide tag, so that test can only ever be indirect.

GRB2 closes part of that gap. MaveDB carries an abundance scan (urn:mavedb:00001305-a-1)
and an activity scan (urn:mavedb:00001305-b-1) built on the SAME library, and that library
contains single-residue deletions and delins variants alongside the substitutions. So the
proxy question can be asked of length-changing variants and of substitutions in one
protein, under one pair of assays, with everything else held constant.

Variant classes, taken from the HGVS protein string:
  missense  p.Ala115Cys
  deletion  p.Ala11del                     (one residue removed)
  delins    p.Ala163_Leu164delinsVal       (two residues replaced by one; net -1)
  nonsense / synonymous: used only as scale anchors, then dropped

Note on the delins class: these are NOT insertions. Two residues are replaced by one, so
the chain gets shorter by one. They are grouped with the deletions as length-changing
variants and reported separately, never described as insertions.

Usage: python 47_indel_proxy.py <mave_dir> <out_json>
"""
import json
import re
import sys

import numpy as np
import pandas as pd
from scipy import stats

RNG = np.random.default_rng(0)
FUNC_DEAD = 0.5
ABUND_OK = 0.75

MIS = re.compile(r"^(?:\w+:)?p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})$")
SYN = re.compile(r"^(?:\w+:)?p\.([A-Z][a-z]{2})(\d+)=$")
DEL = re.compile(r"^(?:\w+:)?p\.([A-Z][a-z]{2})(\d+)del$")
DELINS = re.compile(r"^(?:\w+:)?p\.([A-Z][a-z]{2})(\d+)_([A-Z][a-z]{2})(\d+)delins([A-Z][a-z]{2})$")
TER = re.compile(r"^(?:\w+:)?p\.([A-Z][a-z]{2})(\d+)(?:Ter|\*)$")


def classify(h):
    if not isinstance(h, str):
        return None, None
    h = h.strip()
    for pat, name in ((TER, "nonsense"), (SYN, "synonymous"), (MIS, "missense"),
                      (DEL, "deletion"), (DELINS, "delins")):
        m = pat.match(h)
        if m:
            return name, h
    return None, None


def load(path):
    d = pd.read_csv(path)
    cl = d["hgvs_pro"].map(lambda h: classify(h)[0])
    d = d[cl.notna() & d["score"].notna()].copy()
    d["klass"] = cl[cl.notna()]
    d["key"] = d["hgvs_pro"].str.strip()
    return d


def normalise(d):
    """Use the published WT-relative scale as it stands.

    LABEL-seq publishes both score sets already normalised so that 1 = wild type and
    below 1 = less abundance or activity than wild type, so no rescaling is applied.

    Anchoring on nonsense variants, which is what script 46 does, MUST NOT be used here,
    and the reason is worth recording because it nearly produced a spurious result. In
    this dataset nonsense variants have HIGHER measured abundance (median 1.014) than
    synonymous ones (0.990). The abundance readout therefore does not detect truncation
    at all -- consistent with the abundance signal coming from a terminal tag that a
    premature stop leaves intact. Dividing by (syn - nonsense) = -0.024 inverts and
    explodes the scale, which on a first run turned every correlation in this script
    negative. That was an artefact of the rescaling, not a property of the data.

    The same fact is a real limitation on what this protein can show, and it is carried
    into the report rather than buried: an abundance assay blind to truncation is also
    blind to one whole class of insertion failure.
    """
    d = d.copy()
    d["norm"] = d["score"]
    syn = d.loc[d["klass"] == "synonymous", "score"]
    non = d.loc[d["klass"] == "nonsense", "score"]
    return d, {"scale": "published WT-relative (1 = wild type); no rescaling applied",
               "median_synonymous": round(float(syn.median()), 4),
               "median_nonsense": round(float(non.median()), 4),
               "nonsense_detected_as_loss": bool(non.median() < syn.median()),
               "n_syn": int(len(syn)), "n_nonsense": int(len(non))}


def boot(fn, n, nboot=4000):
    v = [fn(RNG.integers(0, n, n)) for _ in range(nboot)]
    v = [q for q in v if np.isfinite(q)]
    return [round(float(q), 4) for q in np.quantile(v, [.025, .975])] if v else [None, None]


def marks(x, y, label):
    n = len(x)
    rho = float(stats.spearmanr(x, y)[0])
    ab, dead = x >= ABUND_OK, y < FUNC_DEAD
    n_ab = int(ab.sum())
    fp = float((ab & dead).sum() / n_ab) if n_ab else np.nan
    return {"class": label, "n": n,
            "spearman": round(rho, 4),
            "spearman_ci95": boot(lambda i: stats.spearmanr(x[i], y[i])[0], n),
            "n_called_abundant": n_ab,
            "false_pass_rate": round(fp, 4) if np.isfinite(fp) else None,
            "false_pass_ci95": boot(
                lambda i: ((x[i] >= ABUND_OK) & (y[i] < FUNC_DEAD)).sum()
                / max((x[i] >= ABUND_OK).sum(), 1), n),
            "frac_dead": round(float(dead.mean()), 4),
            "median_abundance": round(float(np.median(x)), 4),
            "median_function": round(float(np.median(y)), 4)}


def main(mdir, out_path):
    A, ia = normalise(load(f"{mdir}/GRB2_abundance.csv"))
    F, iff = normalise(load(f"{mdir}/GRB2_function.csv"))
    m = (A[["key", "klass", "norm"]].groupby(["key", "klass"], as_index=False).mean()
         .merge(F[["key", "norm"]].groupby("key", as_index=False).mean(),
                on="key", suffixes=("_ab", "_fn")))

    rep = {"protein": "GRB2",
           "abundance_set": "urn:mavedb:00001305-a-1",
           "function_set": "urn:mavedb:00001305-b-1",
           "scale": {"abundance": ia, "function": iff},
           "note": "delins variants replace two residues with one, so they SHORTEN the "
                   "chain by one; they are length-changing variants, not insertions",
           "by_class": {}}

    for k in ("missense", "deletion", "delins"):
        s = m[m["klass"] == k]
        if len(s) < 30:
            rep["by_class"][k] = {"n": int(len(s)), "status": "too few to score"}
            continue
        rep["by_class"][k] = marks(s["norm_ab"].to_numpy(), s["norm_fn"].to_numpy(), k)

    lc = m[m["klass"].isin(["deletion", "delins"])]
    rep["by_class"]["length_changing_combined"] = marks(
        lc["norm_ab"].to_numpy(), lc["norm_fn"].to_numpy(), "deletion+delins")

    # Usability check, applied after the marks are computed but reported ahead of them,
    # because this dataset fails it and the marks must not be quoted as if it had not.
    fn_all = m["norm_fn"].to_numpy()
    frac_below = float((fn_all < FUNC_DEAD).mean())
    rep["usable_for_this_question"] = {
        "verdict": "NO",
        "reasons": [
            f"only {frac_below:.1%} of variants fall below the loss-of-function threshold "
            f"of {FUNC_DEAD}: the activity assay has almost no dynamic range for loss of "
            "function, so a false-pass rate cannot be estimated from it",
            "nonsense variants score at or above wild type in BOTH assays "
            f"(abundance {ia['median_nonsense']} vs syn {ia['median_synonymous']}; "
            f"activity {iff['median_nonsense']} vs syn {iff['median_synonymous']}), so "
            "neither readout detects a truncated protein as a loss",
            "abundance and activity are ANTI-correlated here (rho about -0.2 in every "
            "class, substitutions included), which is the signature of an overexpression "
            "system in which more GRB2 sequesters signalling components and suppresses "
            "ERK phosphorylation; this makes the pair unusable as a function-vs-abundance "
            "comparison rather than evidence about the proxy question"],
        "consequence": "GRB2 is reported as attempted and rejected. The conclusion rests "
                       "on script 46 (four proteins, substitutions), not on this file."}

    json.dump(rep, open(out_path, "w"), indent=2)
    print("USABLE FOR THE PROXY QUESTION: NO -- see usable_for_this_question in the json")
    print(f"  only {frac_below:.1%} of variants fall below the dead threshold; "
          f"nonsense scores at or above wild type in both assays\n")
    print(f"GRB2  abundance vs function, same library, same assays\n")
    print(f"{'class':24s} {'n':>5s} {'rho':>7s} {'rho 95% CI':>18s} "
          f"{'falsepass':>10s} {'%dead':>6s}")
    for k, v in rep["by_class"].items():
        if "spearman" not in v:
            print(f"{k:24s} {v['n']:5d}  {v.get('status')}")
            continue
        ci = f"[{v['spearman_ci95'][0]:.2f},{v['spearman_ci95'][1]:.2f}]"
        fp = v["false_pass_rate"]
        print(f"{k:24s} {v['n']:5d} {v['spearman']:7.3f} {ci:>18s} "
              f"{(fp if fp is not None else float('nan')):10.3f} {v['frac_dead']:6.3f}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
