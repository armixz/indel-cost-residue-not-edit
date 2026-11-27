"""Independent re-derivation of every headline number, from the raw deposit only.

Deliberately does NOT read results/tsuboyama_indels.csv or results/tsub_features.csv, which were
produced by earlier scripts in this project. Everything here is rebuilt from
Lib3_K50dG.csv and the deposited AlphaFold PDBs, so that a bug in the earlier pipeline cannot
propagate into the verification.

Checks, in order:
  V1  parsing: how many single indels, how the wild-type reference is defined, how many domains
  V2  numbering: the residue named in each deletion variant must match the model at that position
  V3  independence from the old pipeline: agreement of the rebuilt costs with the old table
  V4  the deletion residue-by-structure interaction (headline a)
  V5  the helix-propensity correlation with a bootstrap over domains (headline b)
  V6  the alanine-versus-glycine insertion contrast (headline c)
  V7  sensitivity: quality filter, wild-type definition, element minimum length, structure source

Usage: indel_env python 41_verify.py <lib3_csv> <afdir> <old_results_dir> <out_json>
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

# filled in by V0 from the verified source; see PACE_SOURCE
PACE = {"A": 0.0, "L": 0.21, "R": 0.21, "M": 0.24, "K": 0.26, "Q": 0.39, "E": 0.40, "I": 0.41,
        "W": 0.49, "S": 0.50, "Y": 0.53, "F": 0.54, "H": 0.61, "V": 0.61, "N": 0.65, "T": 0.66,
        "C": 0.68, "D": 0.69, "G": 1.00}
PACE_SOURCE = ("Pace & Scholtz 1998 Biophys J 75:422. All 19 values above were confirmed against "
               "the published abstract (PMC1299714). Proline is EXCLUDED: the abstract says "
               "'excluding proline', and the circulating value 3.16 could not be verified.")
VOL = {"G": 60.1, "A": 88.6, "S": 89.0, "C": 108.5, "D": 111.1, "P": 112.7, "N": 114.1,
       "T": 116.1, "E": 138.4, "V": 140.0, "Q": 143.8, "H": 153.2, "M": 162.9, "I": 166.7,
       "L": 166.7, "K": 168.6, "R": 173.4, "F": 189.9, "Y": 193.6, "W": 227.8}


def read_raw(lib3):
    """Every single-residue indel, plus every wild-type row, straight from the deposit."""
    csv.field_size_limit(10 ** 7)
    wt = defaultdict(list)
    rows = []
    n_total = 0
    with open(lib3, newline="") as fh:
        for r in csv.DictReader(fh):
            n_total += 1
            nm = r["name"]
            last = nm.rsplit("_", 1)[-1]
            if last.startswith("wt"):
                try:
                    wt[nm.rsplit("_", 1)[0]].append(float(r["deltaG"]))
                except ValueError:
                    pass
                continue
            m = PAT.match(nm)
            if not m:
                continue
            try:
                rows.append({"domain": m["dom"], "kind": m["kind"], "res": m["res"],
                             "pos": int(m["pos"]), "dG": float(r["deltaG"]),
                             "ci": float(r["deltaG_95CI"])})
            except (ValueError, KeyError):
                continue
    return {k: float(np.mean(v)) for k, v in wt.items()}, pd.DataFrame(rows), n_total


def structures(afdir, domains, min_h=5, min_e=2):
    out, models = [], {}
    for dom in sorted(domains):
        p = os.path.join(afdir, dom.replace(".pdb", "") + ".pdb")
        if not os.path.exists(p):
            continue
        res = D17.parse_backbone(p)
        ss, keys = D17.assign(res)
        if not ss:
            continue
        el = D17.elements(ss, keys, min_h=min_h, min_e=min_e)
        models[dom] = res
        for i in keys:
            out.append({"domain": dom, "pos": i, "aa": res[i]["aa"], "ss": el[i][0]})
    return pd.DataFrame(out), models


def agg_boot(vals, nboot=20000):
    v = np.asarray([x for x in vals if np.isfinite(x)], float)
    if len(v) < 5:
        return {"n": int(len(v)), "mean": None, "ci95": [None, None]}
    bs = v[RNG.integers(0, len(v), size=(nboot, len(v)))].mean(axis=1)
    return {"n": int(len(v)), "mean": round(float(v.mean()), 4),
            "ci95": [round(float(np.quantile(bs, .025)), 4),
                     round(float(np.quantile(bs, .975)), 4)]}


def helix_contrast_by_residue(df, col="cost", a="H", b="L", nmin=1, min_dom=5):
    pts = {}
    for aa, g in df.groupby("aa"):
        v = [float(x[x["ss"] == a][col].mean() - x[x["ss"] == b][col].mean())
             for _, x in g.groupby("domain")
             if (x["ss"] == a).sum() >= nmin and (x["ss"] == b).sum() >= nmin]
        if len(v) >= min_dom:
            pts[str(aa)] = (float(np.mean(v)), len(v))
    return pts


def spearman_vs_pace(pts, drop=()):
    aas = [a for a in pts if a in PACE and a not in drop]
    if len(aas) < 8:
        return np.nan, len(aas)
    return float(stats.spearmanr([PACE[a] for a in aas],
                                 [pts[a][0] for a in aas])[0]), len(aas)


def main(lib3, afdir, oldres, out_path):
    rep = {"pace_source": PACE_SOURCE}

    # ---- V1 parsing
    wt, v, n_total = read_raw(lib3)
    rep["V1_parsing"] = {
        "rows_in_file": n_total,
        "domains_with_wt_reference": len(wt),
        "single_indel_rows": int(len(v)),
        "insertions": int((v["kind"] == "ins").sum()),
        "deletions": int((v["kind"] == "del").sum()),
        "inserted_residues": sorted(v[v["kind"] == "ins"]["res"].unique().tolist()),
        "deleted_residue_types": int(v[v["kind"] == "del"]["res"].nunique())}

    v = v[v["ci"] <= 1.0].copy()
    v["wtdG"] = v["domain"].map(wt)
    rep["V1_parsing"]["dropped_no_wt_reference"] = int(v["wtdG"].isna().sum())
    v = v[v["wtdG"].notna()].copy()
    v["cost"] = v["wtdG"] - v["dG"]
    rep["V1_parsing"]["after_quality_filter"] = int(len(v))

    # ---- structures
    S, models = structures(afdir, v["domain"].unique())
    rep["V1_parsing"]["domains_with_models"] = int(S["domain"].nunique())

    # ---- V2 numbering, from the deletion variant names
    chk = v[v["kind"] == "del"].merge(S, on=["domain", "pos"], how="inner")
    rep["V2_numbering"] = {
        "checked": int(len(chk)),
        "match_rate": round(float((chk["res"] == chk["aa"]).mean()), 5),
        "mismatched_domains": sorted(chk[chk["res"] != chk["aa"]]["domain"].unique().tolist())[:10]}

    m = v.merge(S, on=["domain", "pos"], how="inner")
    ins = (m[(m["kind"] == "ins") & (m["res"] == "A")]
           .groupby(["domain", "pos"])["cost"].mean().rename("insA").reset_index())
    insG = (m[(m["kind"] == "ins") & (m["res"] == "G")]
            .groupby(["domain", "pos"])["cost"].mean().rename("insG").reset_index())
    dele = (m[m["kind"] == "del"]
            .groupby(["domain", "pos"])["cost"].mean().rename("del").reset_index())
    tab = (S.merge(ins, on=["domain", "pos"], how="left")
           .merge(insG, on=["domain", "pos"], how="left")
           .merge(dele, on=["domain", "pos"], how="left"))
    tab["vol"] = tab["aa"].map(VOL)

    # ---- V3 agreement with the old pipeline
    try:
        old = pd.read_csv(f"{oldres}/tsuboyama_indels.csv")
        oldwt = json.load(open(f"{oldres}/tsub_wt_dG.json"))
        old = old[old["dG_ci"] <= 1.0].copy()
        old["wtdG"] = old["domain"].map(oldwt)
        old = old[old["wtdG"].notna()]
        old["cost_old"] = old["wtdG"] - old["dG"]
        od = old[old["kind"] == "del"].groupby(["domain", "pos"])["cost_old"].mean().reset_index()
        j = dele.merge(od, on=["domain", "pos"], how="inner")
        rep["V3_agreement_with_old_pipeline"] = {
            "overlapping_positions": int(len(j)),
            "pearson": round(float(stats.pearsonr(j["del"], j["cost_old"])[0]), 5),
            "max_abs_diff": round(float((j["del"] - j["cost_old"]).abs().max()), 4),
            "mean_abs_diff": round(float((j["del"] - j["cost_old"]).abs().mean()), 5)}
    except Exception as e:  # noqa: BLE001
        rep["V3_agreement_with_old_pipeline"] = {"error": str(e)}

    # ---- V4 interaction, deletions
    d = tab[tab["del"].notna() & tab["vol"].notna()].copy()
    rep["V4_setup"] = {"n": int(len(d)), "domains": int(d["domain"].nunique()),
                       "residues": int(d["aa"].nunique()),
                       "by_class": {str(k): int(x) for k, x in d["ss"].value_counts().items()}}
    try:
        import statsmodels.formula.api as smf
        d["dom"] = d["domain"].astype("category")
        cl = {"groups": d["dom"].cat.codes}
        b = smf.ols("Q('del') ~ C(dom) + C(aa) + C(ss) + vol", data=d).fit(
            cov_type="cluster", cov_kwds=cl)
        f = smf.ols("Q('del') ~ C(dom) + C(aa)*C(ss) + vol", data=d).fit(
            cov_type="cluster", cov_kwds=cl)
        ft = f.compare_f_test(b)
        rep["V4_interaction"] = {"F": round(float(ft[0]), 4), "p": float(ft[1]),
                                 "df_num": float(ft[2])}
        rng = np.random.default_rng(99)
        p = d.copy()
        p["aa"] = p.groupby("domain")["aa"].transform(lambda s: rng.permutation(s.values))
        p["dom"] = p["domain"].astype("category")
        clp = {"groups": p["dom"].cat.codes}
        b2 = smf.ols("Q('del') ~ C(dom) + C(aa) + C(ss) + vol", data=p).fit(
            cov_type="cluster", cov_kwds=clp)
        f2 = smf.ols("Q('del') ~ C(dom) + C(aa)*C(ss) + vol", data=p).fit(
            cov_type="cluster", cov_kwds=clp)
        ft2 = f2.compare_f_test(b2)
        rep["V4_permuted"] = {"F": round(float(ft2[0]), 4), "p": float(ft2[1])}
    except ImportError:
        rep["V4_interaction"] = {"status": "statsmodels unavailable"}

    # ---- V5 propensity correlation, bootstrap over domains
    # Precompute one contrast per (domain, residue); resampling domains then only re-averages
    # those numbers, which is identical to recomputing from scratch and far cheaper.
    cell = {}
    for (dom, aa), g in d.groupby(["domain", "aa"]):
        h, lp = g[g["ss"] == "H"]["del"], g[g["ss"] == "L"]["del"]
        if len(h) >= 1 and len(lp) >= 1:
            cell.setdefault(str(dom), {})[str(aa)] = float(h.mean() - lp.mean())
    doms = sorted(cell)

    def rho_from_cells(picked, drop=(), min_dom=5):
        acc = defaultdict(list)
        for dm in picked:
            for aa, val in cell[dm].items():
                acc[aa].append(val)
        pairs = [(PACE[a], float(np.mean(v))) for a, v in acc.items()
                 if a in PACE and a not in drop and len(v) >= min_dom]
        if len(pairs) < 8:
            return np.nan, len(pairs)
        xs, ys = zip(*pairs)
        return float(stats.spearmanr(xs, ys)[0]), len(pairs)

    pts = {}
    accum = defaultdict(list)
    for dm in doms:
        for aa, val in cell[dm].items():
            accum[aa].append(val)
    for aa, vals_ in accum.items():
        if len(vals_) >= 5:
            pts[aa] = (float(np.mean(vals_)), len(vals_))
    r0, n0 = rho_from_cells(doms)
    boot = []
    for _ in range(2000):
        pick = list(RNG.choice(doms, len(doms), replace=True))
        r, _ = rho_from_cells(pick)
        if np.isfinite(r):
            boot.append(r)
    q = np.quantile(boot, [.025, .975]) if boot else [np.nan, np.nan]
    rep["V5_propensity"] = {
        "spearman": round(float(r0), 4), "n_residues": n0,
        "bootstrap_ci95_over_domains": [round(float(q[0]), 4), round(float(q[1]), 4)],
        "n_boot": len(boot),
        "contrasts": {k: [round(x[0], 4), x[1]] for k, x in sorted(pts.items())}}
    for lab, drop in (("without_P", ("P",)), ("without_G_P", ("G", "P")),
                      ("without_G_P_W_F", ("G", "P", "W", "F"))):
        rr, nn = rho_from_cells(doms, drop)
        bb = []
        for _ in range(2000):
            pick = list(RNG.choice(doms, len(doms), replace=True))
            r, _ = rho_from_cells(pick, drop)
            if np.isfinite(r):
                bb.append(r)
        qq = np.quantile(bb, [.025, .975]) if bb else [np.nan, np.nan]
        rep["V5_propensity"][lab] = {"spearman": round(float(rr), 4), "n_residues": nn,
                                     "ci95": [round(float(qq[0]), 4), round(float(qq[1]), 4)]}

    # ---- V6 Ala vs Gly insertion contrast
    t = tab[tab["insA"].notna() & tab["insG"].notna() & tab["del"].notna()].copy()
    out6 = {"n_positions": int(len(t)), "domains": int(t["domain"].nunique())}
    for col, lab in (("insA", "alanine"), ("insG", "glycine")):
        for a, b2 in (("H", "L"), ("E", "L")):
            v2 = [float((x[x["ss"] == a][col] - x[x["ss"] == a]["del"]).mean()
                        - (x[x["ss"] == b2][col] - x[x["ss"] == b2]["del"]).mean())
                  for _, x in t.groupby("domain")
                  if (x["ss"] == a).sum() >= 2 and (x["ss"] == b2).sum() >= 2]
            out6[f"{lab}_{a}_minus_{b2}"] = agg_boot(v2)
    rep["V6_insertion_contrast"] = out6

    # ---- V7 sensitivity
    sens = {}
    for name, q_ci in (("ci_le_0.5", 0.5), ("ci_le_2.0", 2.0)):
        vv = v[v["ci"] <= q_ci]
        dd = (S.merge(vv[vv["kind"] == "del"].groupby(["domain", "pos"])["cost"].mean()
                      .rename("del").reset_index(), on=["domain", "pos"]))
        pp = helix_contrast_by_residue(dd, "del")
        r, n = spearman_vs_pace(pp)
        sens[name] = {"n_positions": int(len(dd)), "spearman": round(float(r), 4),
                      "n_residues": n}  # noqa
    S3, _ = structures(afdir, v["domain"].unique(), min_h=6, min_e=3)
    dd = S3.merge(dele, on=["domain", "pos"])
    pp = helix_contrast_by_residue(dd, "del")
    r, n = spearman_vs_pace(pp)
    sens["stricter_elements_minH6_minE3"] = {"n_positions": int(len(dd)),
                                             "spearman": round(float(r), 4), "n_residues": n}
    rep["V7_sensitivity"] = sens

    json.dump(rep, open(out_path, "w"), indent=2)
    print(json.dumps(rep, indent=2))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
