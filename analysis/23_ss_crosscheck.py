"""Cross-check 17_dssp.py against Topolska et al.'s own secondary-structure table.

`tsuboyama_sec_struc.rds` in Zenodo 10830272 (additional_files) holds the Lehner lab's structure
assignment for the Tsuboyama domains: 13,682 positions with Structure, simple_struc, element
position, element length and rSASA. It was produced independently of this project, so agreement
with 17_dssp.py is external evidence that the assignment used here is sound, and disagreement is
a warning.

Checks:
  1. three-state agreement (H / E / L) position by position
  2. agreement of element length for positions both call helix or strand
  3. correlation of their rSASA with the rsa computed here
  4. whether the key result -- the insertion/deletion asymmetry by structure class -- is the same
     when THEIR labels are used instead of ours

Usage: indel_env python 23_ss_crosscheck.py <results_dir> <rds_path> <out_json>
"""
import json
import sys

import numpy as np
import pandas as pd
import pyreadr
from scipy import stats

CI_MAX = 1.0
RNG = np.random.default_rng(0)

MAP = {"AlphaHelix": "H", "Helix": "H", "310Helix": "H", "PiHelix": "H",
       "Strand": "E", "Bridge": "E", "BetaBridge": "E",
       "loop": "L", "Loop": "L", "Turn": "L", "Bend": "L",
       "ntermini": "T", "ctermini": "T"}


def agg(vals, nboot=10000):
    v = np.asarray([x for x in vals if x is not None and np.isfinite(x)], dtype=float)
    if len(v) < 10:
        return {"n_domains": int(len(v)), "mean": None, "ci95": [None, None]}
    idx = RNG.integers(0, len(v), size=(nboot, len(v)))
    bs = v[idx].mean(axis=1)
    return {"n_domains": int(len(v)), "mean": round(float(v.mean()), 3),
            "ci95": [round(float(np.quantile(bs, 0.025)), 3),
                     round(float(np.quantile(bs, 0.975)), 3)]}


def main(rdir, rds, out_path):
    their = list(pyreadr.read_r(rds).values())[0]
    their = their.rename(columns={"Pos": "struct_pos", "pdb_name": "domain"})
    their["struct_pos"] = their["struct_pos"].astype(int)
    their["ss_them"] = their["simple_struc"].map(MAP).fillna("L")
    mine = pd.read_csv(f"{rdir}/tsub_ss.csv")
    mine["domain"] = mine["domain"].str.replace(".pdb", "", regex=False)

    m = mine.merge(their[["domain", "struct_pos", "resid", "ss_them", "simple_struc",
                          "element_lenght", "rSASA"]],
                   on=["domain", "struct_pos"], how="inner")
    rep = {"merged_positions": int(len(m)),
           "domains_in_common": int(m["domain"].nunique()),
           "residue_identity_match": round(float((m["aa"] == m["resid"]).mean()), 4)}

    # termini are not a structure class here; compare only where they call H/E/L
    core = m[m["ss_them"].isin(["H", "E", "L"])]
    rep["three_state_agreement"] = round(float((core["ss"] == core["ss_them"]).mean()), 4)
    rep["confusion"] = json.loads(
        pd.crosstab(core["ss_them"], core["ss"]).to_json(orient="index"))
    rep["their_fractions"] = {k: round(float(v), 4) for k, v in
                              core["ss_them"].value_counts(normalize=True).items()}
    rep["my_fractions"] = {k: round(float(v), 4) for k, v in
                           core["ss"].value_counts(normalize=True).items()}

    both = core[(core["ss"] == core["ss_them"]) & (core["ss"].isin(["H", "E"]))]
    rep["element_length_agreement"] = {
        "n": int(len(both)),
        "pearson": round(float(stats.pearsonr(both["elem_len"], both["element_lenght"])[0]), 3),
        "median_abs_diff": round(float((both["elem_len"] - both["element_lenght"]).abs().median()), 2)}

    r = m[m["rSASA"].notna()]
    feat = pd.read_csv(f"{rdir}/tsub_features.csv")
    r = r.merge(feat[["domain", "struct_pos", "rsa"]].drop_duplicates(),
                on=["domain", "struct_pos"], how="inner")
    rep["rsa_agreement"] = {"n": int(len(r)),
                            "pearson": round(float(stats.pearsonr(r["rsa"], r["rSASA"])[0]), 3)}

    # does the headline asymmetry survive with THEIR labels?
    ind = pd.read_csv(f"{rdir}/tsuboyama_indels.csv")
    wt = json.load(open(f"{rdir}/tsub_wt_dG.json"))
    ind = ind[ind["dG_ci"] <= CI_MAX].copy()
    ind["wt"] = ind["domain"].map(wt)
    ind = ind[ind["wt"].notna()]
    ind["cost"] = ind["wt"] - ind["dG"]
    tab = feat[["domain", "pos", "struct_pos"]].copy()
    for key, kind, res in (("ins_A", "ins", "A"), ("del", "del", None)):
        sub = ind[ind["kind"] == kind]
        if res is not None:
            sub = sub[sub["res"] == res]
        t = sub.groupby(["domain", "pos"])["cost"].mean().rename(key).reset_index()
        tab = tab.merge(t, on=["domain", "pos"], how="left")
    tab = tab[tab["ins_A"].notna() & tab["del"].notna()].copy()
    tab["asym"] = tab["ins_A"] - tab["del"]
    j = tab.merge(core[["domain", "struct_pos", "ss", "ss_them"]], on=["domain", "struct_pos"])
    rep["asymmetry_by_their_labels"] = {
        s: agg([float(x["asym"].mean()) for _, x in g.groupby("domain")])
        for s, g in j.groupby("ss_them")}
    rep["asymmetry_by_my_labels_same_rows"] = {
        s: agg([float(x["asym"].mean()) for _, x in g.groupby("domain")])
        for s, g in j.groupby("ss")}

    json.dump(rep, open(out_path, "w"), indent=2)
    print(json.dumps(rep, indent=2))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])
