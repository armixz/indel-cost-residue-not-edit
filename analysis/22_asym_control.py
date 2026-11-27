"""Is the insertion-vs-deletion asymmetry a residue-identity artefact?

H3 of PREREG_indel.md (declared exploratory) found that insertions are cheaper than deletions
everywhere, and twice as much in beta strands (-0.587) as in helices (-0.266) or loops (-0.243).

A deletion always removes the WILD-TYPE residue, while the insertion is always Ala or Gly. Strands
are rich in large buried hydrophobics, so the strand excess could be nothing but "deleting a large
buried residue is expensive". This script tries to kill that explanation before any mechanistic
story is entertained. It is a control, not a new hypothesis.

Arms, in increasing strength:
  1. raw asymmetry by secondary structure (reproduces H3)
  2. asymmetry restricted to positions whose wild-type residue IS alanine (identity matched to the
     Ala insertion, so "what" is held fixed and only "where" varies)
  3. asymmetry stratified by wild-type residue volume tertile
  4. asymmetry stratified by burial (RSA tertile), within strands and helices separately
  5. per-residue-type asymmetry, to see whether the strand excess survives inside residue classes
  6. a direct regression check: does ss still carry the asymmetry once volume and RSA are in?

All statistics are computed within domain and then averaged over domains, as in 21_indel_marks2.py.

Usage: run.sh 22_asym_control.py <results_dir> <out_json>
"""
import json
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, "/mnt/c/Users/armin/Desktop/AXZ/insertion-lab/analysis")

CI_MAX = 1.0
NBOOT = 10000
RNG = np.random.default_rng(0)

# Zamyatnin residue volumes, cubic angstrom
VOL = {"G": 60.1, "A": 88.6, "S": 89.0, "C": 108.5, "D": 111.1, "P": 112.7, "N": 114.1,
       "T": 116.1, "E": 138.4, "V": 140.0, "Q": 143.8, "H": 153.2, "M": 162.9, "I": 166.7,
       "L": 166.7, "K": 168.6, "R": 173.4, "F": 189.9, "Y": 193.6, "W": 227.8}


def agg(vals, nboot=NBOOT):
    v = np.asarray([x for x in vals if x is not None and np.isfinite(x)], dtype=float)
    if len(v) < 10:
        return {"n_domains": int(len(v)), "mean": None, "ci95": [None, None]}
    idx = RNG.integers(0, len(v), size=(nboot, len(v)))
    bs = v[idx].mean(axis=1)
    return {"n_domains": int(len(v)), "n_used": int(len(v)), "mean": round(float(v.mean()), 3),
            "ci95": [round(float(np.quantile(bs, 0.025)), 3),
                     round(float(np.quantile(bs, 0.975)), 3)]}


def dom_mean(df, col="asym"):
    return agg([float(g[col].mean()) for _, g in df.groupby("domain") if len(g)])


def load(rdir):
    ind = pd.read_csv(f"{rdir}/tsuboyama_indels.csv")
    feat = pd.read_csv(f"{rdir}/tsub_features.csv")
    ss = pd.read_csv(f"{rdir}/tsub_ss.csv")
    ss["domain"] = ss["domain"].str.replace(".pdb", "", regex=False)
    wt = json.load(open(f"{rdir}/tsub_wt_dG.json"))
    ind = ind[ind["dG_ci"] <= CI_MAX].copy()
    ind["wt"] = ind["domain"].map(wt)
    ind = ind[ind["wt"].notna()]
    ind["cost"] = ind["wt"] - ind["dG"]
    df = feat[["domain", "pos", "struct_pos", "plddt", "rsa"]].merge(
        ss[["domain", "struct_pos", "aa", "ss", "idx_end", "elem_len", "d_loop"]],
        on=["domain", "struct_pos"], how="inner")
    for key, kind, res in (("ins_A", "ins", "A"), ("ins_G", "ins", "G"), ("del", "del", None)):
        sub = ind[ind["kind"] == kind]
        if res is not None:
            sub = sub[sub["res"] == res]
        t = sub.groupby(["domain", "pos"])["cost"].mean().rename(key).reset_index()
        df = df.merge(t, on=["domain", "pos"], how="left")
    df = df[df["ins_A"].notna() & df["del"].notna()].copy()
    df["asym"] = df["ins_A"] - df["del"]
    df["asym_G"] = df["ins_G"] - df["del"]
    df["vol"] = df["aa"].map(VOL)
    return df


def main(rdir, out_path):
    df = load(rdir)
    rep = {"n_positions": int(len(df)), "n_domains": int(df["domain"].nunique())}

    # 1. raw
    rep["arm1_raw_by_ss"] = {str(s): dom_mean(g) for s, g in df.groupby("ss")}

    # 2. wild-type alanine only -- identity matched
    ala = df[df["aa"] == "A"]
    rep["arm2_wt_alanine_only"] = {
        "n_positions": int(len(ala)),
        "by_ss": {str(s): dom_mean(g) for s, g in ala.groupby("ss")},
        "note": ("the insertion is Ala and the deleted residue is Ala, so chemistry is matched; "
                 "any asymmetry left is about where, not what")}

    # 3. volume tertiles
    d = df[df["vol"].notna()].copy()
    d["vt"] = pd.qcut(d["vol"], 3, labels=["small", "mid", "large"], duplicates="drop")
    rep["arm3_by_volume"] = {str(t): {"mean_vol": round(float(g["vol"].mean()), 1),
                                      "all": dom_mean(g),
                                      "by_ss": {str(s): dom_mean(x) for s, x in g.groupby("ss")}}
                             for t, g in d.groupby("vt", observed=True)}

    # 4. burial tertiles within E and H
    r = df[df["rsa"].notna()].copy()
    r["bt"] = pd.qcut(r["rsa"], 3, labels=["buried", "mid", "exposed"], duplicates="drop")
    rep["arm4_by_burial"] = {}
    for s in ("E", "H", "L"):
        sub = r[r["ss"] == s]
        rep["arm4_by_burial"][s] = {str(t): dom_mean(g)
                                    for t, g in sub.groupby("bt", observed=True)}

    # 5. per wild-type residue, strand minus the rest
    per_aa = {}
    for a, g in df.groupby("aa"):
        if len(g) < 100:
            continue
        e = g[g["ss"] == "E"]
        o = g[g["ss"] != "E"]
        if len(e) < 30 or len(o) < 30:
            continue
        per_aa[str(a)] = {"n_strand": int(len(e)), "strand": dom_mean(e), "other": dom_mean(o),
                          "strand_minus_other": round(float(e["asym"].mean() - o["asym"].mean()), 3)}
    rep["arm5_per_wt_residue"] = per_aa
    diffs = [v["strand_minus_other"] for v in per_aa.values()]
    rep["arm5_summary"] = {
        "n_residue_types": len(diffs),
        "n_negative": int(sum(1 for x in diffs if x < 0)),
        "median_strand_minus_other": round(float(np.median(diffs)), 3) if diffs else None,
        "note": "if the strand excess is only about which residue is deleted, it should vanish "
                "inside a fixed residue type"}

    # 6. regression: does ss survive volume + rsa + plddt?
    try:
        import statsmodels.api as sm
        x = df[df["vol"].notna() & df["rsa"].notna() & df["plddt"].notna()].copy()
        X = pd.DataFrame({"vol": x["vol"], "rsa": x["rsa"], "plddt": x["plddt"],
                          "is_E": (x["ss"] == "E").astype(float),
                          "is_H": (x["ss"] == "H").astype(float)})
        X = sm.add_constant(X)
        m = sm.OLS(x["asym"], X).fit(cov_type="cluster",
                                     cov_kwds={"groups": x["domain"].astype("category").cat.codes})
        rep["arm6_regression"] = {"n": int(len(x)),
                                  "coef": {k: round(float(v), 4) for k, v in m.params.items()},
                                  "ci95": {k: [round(float(a), 4), round(float(b), 4)]
                                           for k, (a, b) in m.conf_int().iterrows()},
                                  "note": "OLS with domain-clustered standard errors; is_E is the "
                                          "strand excess after volume, burial and pLDDT"}
    except ImportError:
        rep["arm6_regression"] = {"status": "statsmodels not available"}

    # Gly arm of the same control
    dfg = df[df["ins_G"].notna()].copy()
    rep["gly_arm_raw_by_ss"] = {str(s): dom_mean(g, "asym_G") for s, g in dfg.groupby("ss")}
    alag = dfg[dfg["aa"] == "G"]
    rep["gly_arm_wt_glycine_only"] = {"n_positions": int(len(alag)),
                                      "by_ss": {str(s): dom_mean(g, "asym_G")
                                                for s, g in alag.groupby("ss")}}

    json.dump(rep, open(out_path, "w"), indent=2)
    print(json.dumps(rep, indent=2))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
