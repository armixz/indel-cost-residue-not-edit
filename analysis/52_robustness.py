"""The full robustness battery, pre-registered in PREREG_robustness.md.

Run once, as one batch, so that the manuscript can be edited once and frozen. Every mark
defined in the pre-registration is computed here and written to the output, whether it
supports the manuscript or not.

R0  composition of the analysis set (a correction, not a test)
R1  natural versus designed replication
R2  a second, independent helix propensity scale (O'Neil & DeGrado 1990)
R3  the variant quality filter sweep
R4  estimator choices: minimum-domain threshold, mean versus median
R5  influence: leave-one-residue-out, leave-one-domain-out, glycine excluded
R6  a permutation null for the propensity correlation itself

Usage: python 52_robustness.py <lib3_csv> <afdir> <results_dir> <out_json>
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
sp = importlib.util.spec_from_file_location("d43", os.path.join(HERE, "43_discriminate.py"))
D43 = importlib.util.module_from_spec(sp)
sp.loader.exec_module(D43)

RNG = np.random.default_rng(0)
PACE = dict(D43.PACE)

# O'Neil & DeGrado, Science 1990;250:646, Table 2, column ddG_alpha (kcal/mol).
# Read by me from a photographic reproduction of the printed table. As printed the
# convention is already "larger = worse helix former": Ala -0.77 is the best, Gly 0.00 is
# the reference, Pro is about +3. Aib is omitted as non-standard. Proline is omitted here
# for the same reason it is omitted from the Pace analysis, and its listed value of about
# 3 is the likely origin of the unverifiable "3.16" that circulates for the Pace scale.
ONEIL = {"A": -0.77, "R": -0.68, "K": -0.65, "L": -0.62, "M": -0.50, "W": -0.45,
         "F": -0.41, "S": -0.35, "Q": -0.33, "E": -0.27, "C": -0.23, "I": -0.23,
         "Y": -0.17, "D": -0.15, "V": -0.14, "T": -0.11, "N": -0.07, "H": -0.06,
         "G": 0.00}

NAT = re.compile(r"^[0-9][A-Za-z0-9]{3}\.pdb$")
DES = re.compile(r"(HEEH|EEHEE|HHH|EHEE|_rd[0-9]|TrROS|Hall)")


def kind(dom):
    if NAT.match(dom):
        return "natural"
    if DES.search(dom):
        return "designed"
    return "other"


def load_raw(lib3, afdir, rdir, ci_max):
    """D43.load with a configurable quality filter. Mirrors its body exactly otherwise."""
    csv.field_size_limit(10 ** 7)
    wt, rows = {}, []
    with open(lib3, newline="") as fh:
        for r in csv.DictReader(fh):
            nm = r["name"]
            if nm.rsplit("_", 1)[-1].startswith("wt"):
                wt.setdefault(nm.rsplit("_", 1)[0], []).append(float(r["deltaG"]))
                continue
            m = D43.PAT.match(nm)
            if not m or m["kind"] != "del":
                continue
            try:
                rows.append({"domain": m["dom"], "aa": m["res"], "pos": int(m["pos"]),
                             "dG": float(r["deltaG"]), "ci": float(r["deltaG_95CI"])})
            except (ValueError, KeyError):
                continue
    wtm = {k: float(np.mean(v)) for k, v in wt.items()}
    v = pd.DataFrame(rows)
    if ci_max is not None:
        v = v[v["ci"] <= ci_max]
    v = v.copy()
    v["wtdG"] = v["domain"].map(wtm)
    v = v[v["wtdG"].notna()]
    v["cost"] = v["wtdG"] - v["dG"]
    S = pd.read_csv(os.path.join(rdir, "tsub_ss.csv"))[
        ["domain", "struct_pos", "aa", "ss"]].rename(
        columns={"struct_pos": "pos", "aa": "aa_model"})
    d = v.groupby(["domain", "pos", "aa"])["cost"].mean().reset_index().merge(
        S, on=["domain", "pos"])
    d = d[d["aa"] == d["aa_model"]]
    return d[d["aa"].isin(PACE)].copy()


def contrasts(d, a="H", b="L", min_dom=5, agg="mean"):
    cell = {}
    for (dom, aa), g in d.groupby(["domain", "aa"]):
        x, y = g[g["ss"] == a]["cost"], g[g["ss"] == b]["cost"]
        if len(x) >= 1 and len(y) >= 1:
            cell.setdefault(str(aa), []).append(float(x.mean() - y.mean()))
    f = np.mean if agg == "mean" else np.median
    return {aa: float(f(v)) for aa, v in cell.items() if len(v) >= min_dom}


def rho_of(ctr, scale=None):
    scale = scale or PACE
    k = sorted(a for a in ctr if a in scale)
    if len(k) < 8:
        return np.nan, 0
    return float(stats.spearmanr([scale[a] for a in k], [ctr[a] for a in k])[0]), len(k)


def boot_rho(d, scale=None, nboot=1000, min_dom=5):
    doms = sorted(d["domain"].unique())
    cell = {}
    for (dom, aa), g in d.groupby(["domain", "aa"]):
        x, y = g[g["ss"] == "H"]["cost"], g[g["ss"] == "L"]["cost"]
        if len(x) >= 1 and len(y) >= 1:
            cell.setdefault(str(dom), {})[str(aa)] = float(x.mean() - y.mean())
    vals = []
    for _ in range(nboot):
        pick = RNG.choice(doms, len(doms), replace=True)
        acc = {}
        for dm in pick:
            for aa, v in cell.get(dm, {}).items():
                acc.setdefault(aa, []).append(v)
        ctr = {a: float(np.mean(v)) for a, v in acc.items() if len(v) >= min_dom}
        r, _ = rho_of(ctr, scale)
        if np.isfinite(r):
            vals.append(r)
    if not vals:
        return [None, None]
    return [round(float(q), 4) for q in np.quantile(vals, [.025, .975])]


def main(lib3, afdir, rdir, out_path):
    rep = {"prereg": "PREREG_robustness.md"}
    d = load_raw(lib3, afdir, rdir, 1.0)
    d["kind"] = d["domain"].map(kind)

    # ---- R0 composition
    comp = d.groupby("kind")["domain"].nunique().to_dict()
    rep["R0_composition"] = {
        "domains_total": int(d["domain"].nunique()),
        "domains_by_kind": {k: int(v) for k, v in comp.items()},
        "rows_by_kind": {k: int(v) for k, v in d.groupby("kind").size().items()},
        "note": "the manuscript described this set as 748 natural domains; it is not"}

    base_ctr = contrasts(d)
    base_rho, base_n = rho_of(base_ctr)
    rep["baseline"] = {"spearman": round(base_rho, 4), "n_residues": base_n,
                       "ci95": boot_rho(d)}

    # ---- R1 natural vs designed
    rep["R1_natural_vs_designed"] = {}
    for k in ("natural", "designed", "other"):
        sub = d[d["kind"] == k]
        if sub["domain"].nunique() < 20:
            rep["R1_natural_vs_designed"][k] = {"domains": int(sub["domain"].nunique()),
                                                "status": "too few"}
            continue
        c = contrasts(sub)
        r, n = rho_of(c)
        rep["R1_natural_vs_designed"][k] = {
            "domains": int(sub["domain"].nunique()), "rows": int(len(sub)),
            "n_residues": n, "spearman": round(r, 4), "ci95": boot_rho(sub),
            "mean_contrast": round(float(np.mean(list(c.values()))), 4)}

    # R1b interaction F test per subset
    try:
        import statsmodels.formula.api as smf
        rep["R1b_interaction_F"] = {}
        for k in ("natural", "designed"):
            sub = d[d["kind"] == k].copy()
            sub["pace"] = sub["aa"].map(PACE)
            sub["isH"] = (sub["ss"] == "H").astype(float)
            sub["isE"] = (sub["ss"] == "E").astype(float)
            m0 = smf.ols("cost ~ C(domain) + C(aa) + isH + isE", data=sub).fit()
            m1 = smf.ols("cost ~ C(domain) + C(aa)*isH + C(aa)*isE", data=sub).fit()
            cmp_ = m1.compare_f_test(m0)
            rep["R1b_interaction_F"][k] = {
                "F": round(float(cmp_[0]), 3), "p": float(cmp_[1]),
                "n": int(len(sub)), "domains": int(sub["domain"].nunique())}
    except Exception as e:  # noqa: BLE001
        rep["R1b_interaction_F"] = {"status": f"unavailable: {e}"}

    # ---- R2 second helix scale
    r2, n2 = rho_of(base_ctr, ONEIL)
    k = sorted(set(PACE) & set(ONEIL))
    rep["R2_second_scale"] = {
        "oneil_degrado": {"spearman": round(r2, 4), "n_residues": n2,
                          "ci95": boot_rho(d, ONEIL)},
        "pace_vs_oneil_scale_correlation": round(float(stats.spearmanr(
            [PACE[a] for a in k], [ONEIL[a] for a in k])[0]), 4),
        "note": "O'Neil & DeGrado list proline at about 3 kcal/mol, the likely origin of "
                "the unverifiable 3.16 that circulates for the Pace scale"}

    # ---- R3 quality filter sweep
    rep["R3_quality_filter"] = {}
    for thr in (0.5, 0.75, 1.0, 1.5, 2.0, None):
        dd = load_raw(lib3, afdir, rdir, thr)
        c = contrasts(dd)
        r, n = rho_of(c)
        rep["R3_quality_filter"][str(thr)] = {
            "rows": int(len(dd)), "domains": int(dd["domain"].nunique()),
            "n_residues": n, "spearman": round(r, 4) if np.isfinite(r) else None}

    # ---- R4 estimator
    rep["R4_estimator"] = {"min_domains": {}, "aggregator": {}}
    for md in (3, 5, 8, 10):
        c = contrasts(d, min_dom=md)
        r, n = rho_of(c)
        rep["R4_estimator"]["min_domains"][str(md)] = {
            "n_residues": n, "spearman": round(r, 4) if np.isfinite(r) else None}
    for agg in ("mean", "median"):
        c = contrasts(d, agg=agg)
        r, n = rho_of(c)
        rep["R4_estimator"]["aggregator"][agg] = {
            "n_residues": n, "spearman": round(r, 4)}

    # ---- R5 influence
    loo_res = {}
    for aa in sorted(base_ctr):
        c = {k_: v for k_, v in base_ctr.items() if k_ != aa}
        r, _ = rho_of(c)
        loo_res[aa] = round(r, 4)
    worst = max(loo_res, key=lambda a: abs(loo_res[a] - base_rho))
    doms = sorted(d["domain"].unique())
    loo_dom = []
    for dm in doms:
        c = contrasts(d[d["domain"] != dm])
        r, _ = rho_of(c)
        if np.isfinite(r):
            loo_dom.append(r)
    c_nog = {k_: v for k_, v in base_ctr.items() if k_ != "G"}
    r_nog, n_nog = rho_of(c_nog)
    rep["R5_influence"] = {
        "leave_one_residue_out": loo_res,
        "most_influential_residue": {"residue": worst, "spearman_without_it": loo_res[worst],
                                     "shift": round(loo_res[worst] - base_rho, 4)},
        "leave_one_domain_out": {"n_domains": len(loo_dom),
                                 "min": round(float(min(loo_dom)), 4),
                                 "max": round(float(max(loo_dom)), 4)},
        "without_glycine": {"spearman": round(r_nog, 4), "n_residues": n_nog}}

    # ---- R6 permutation null
    null = []
    for _ in range(500):
        p = d.copy()
        p["aa"] = p.groupby("domain")["aa"].transform(
            lambda s: RNG.permutation(s.to_numpy()))
        c = contrasts(p)
        r, _ = rho_of(c)
        if np.isfinite(r):
            null.append(r)
    rep["R6_permutation_null"] = {
        "n_permutations": len(null),
        "mean": round(float(np.mean(null)), 4),
        "sd": round(float(np.std(null)), 4),
        "quantiles_2.5_97.5": [round(float(q), 4) for q in np.quantile(null, [.025, .975])],
        "observed": round(base_rho, 4),
        "observed_more_extreme_than_all_permutations":
            bool(base_rho < min(null)) if null else None}

    json.dump(rep, open(out_path, "w"), indent=2)
    print(json.dumps(rep, indent=2)[:400])
    print("\n... full report written to", out_path)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
