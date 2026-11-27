"""Re-computation of Gate 11 with hydrogen-bond secondary structure.

NOT a new hypothesis. Gate 11 (RESULTS.md) reported ddG(Ala) - ddG(Gly) of +0.201 in helix,
+0.102 in strand and -0.082 in loop, helix minus loop +0.283 with by-domain CI [0.244, 0.321].
Those classes came from 16_gate11_gly_ala.py's ss3(), which labels any residue inside two very
wide Ramachandran boxes as H or E and leaves only ~10% of residues as loop. 17_dssp.py replaces
it with Kabsch-Sander hydrogen-bond assignment, validated in 18_ss_check.py against dihedral
geometry it never used (H: phi -63.0 psi -41.8, CA(i,i+4) 6.2 A).

This script recomputes the SAME statistic as Gate 11 with the corrected classes, so the two are
directly comparable. The marks are Gate 11's own, unchanged:
  A: helix minus loop > 0, by-domain CI excluding 0
  B: helix minus loop >= 0.3 kcal/mol

Usage: run.sh 19_gate11_dssp.py <results_dir> <out_json>
"""
import json
import sys

import numpy as np
import pandas as pd

CI_MAX = 1.0  # same filter Gate 11 used


def boot_by_domain(df, fn, n=2000, seed=0):
    rng = np.random.default_rng(seed)
    doms = df["domain"].unique()
    idx = {d: df.index[df["domain"] == d].to_numpy() for d in doms}
    vals = []
    for _ in range(n):
        pick = rng.choice(doms, len(doms), replace=True)
        sub = df.loc[np.concatenate([idx[d] for d in pick])]
        v = fn(sub)
        if v is not None and np.isfinite(v):
            vals.append(v)
    return [round(float(q), 3) for q in np.nanquantile(vals, [0.025, 0.975])]


def main(rdir, out_path):
    ind = pd.read_csv(f"{rdir}/tsuboyama_indels.csv")
    feat = pd.read_csv(f"{rdir}/tsub_features.csv")
    ss = pd.read_csv(f"{rdir}/tsub_ss.csv")
    ss["domain"] = ss["domain"].str.replace(".pdb", "", regex=False)
    wt = json.load(open(f"{rdir}/tsub_wt_dG.json"))

    ins = ind[(ind["kind"] == "ins") & (ind["dG_ci"] <= CI_MAX)].copy()
    ins["wt"] = ins["domain"].map(wt)
    ins = ins[ins["wt"].notna()]
    ins["ddg"] = ins["dG"] - ins["wt"]
    piv = (ins.pivot_table(index=["domain", "pos"], columns="res", values="ddg")
           .dropna().reset_index())
    piv["delta"] = piv["A"] - piv["G"]

    df = (piv.merge(feat[["domain", "pos", "struct_pos", "plddt", "rsa"]], on=["domain", "pos"])
          .merge(ss[["domain", "struct_pos", "ss", "idx_end", "elem_len", "d_loop"]],
                 on=["domain", "struct_pos"]))

    rep = {"n_positions": int(len(df)), "n_domains": int(df["domain"].nunique())}
    by = df.groupby("ss")["delta"].agg(["count", "mean", "median"]).round(3)
    rep["by_ss"] = json.loads(by.to_json(orient="index"))

    hl = lambda s: (s[s["ss"] == "H"]["delta"].mean() - s[s["ss"] == "L"]["delta"].mean())
    diff = hl(df)
    ci = boot_by_domain(df, hl)
    rep["A_helix_minus_loop"] = {"difference_kcal": round(float(diff), 3), "ci95_by_domain": ci,
                                 "mark_A_pass": bool(ci[0] > 0),
                                 "mark_B_pass": bool(diff >= 0.3)}

    el = lambda s: (s[s["ss"] == "E"]["delta"].mean() - s[s["ss"] == "L"]["delta"].mean())
    rep["strand_minus_loop"] = {"difference_kcal": round(float(el(df)), 3),
                                "ci95_by_domain": boot_by_domain(df, el)}

    # helix middle (>=4 from both ends) against helix ends (idx_end <= 1), post-hoc in Gate 11
    h = df[df["ss"] == "H"].copy()
    h["where"] = np.where(h["idx_end"] >= 4, "middle", np.where(h["idx_end"] <= 1, "end", "mid1"))
    rep["helix_middle_vs_end"] = json.loads(
        h.groupby("where")["delta"].agg(["count", "mean"]).round(3).to_json(orient="index"))
    mv = lambda s: (s[(s["ss"] == "H") & (s["idx_end"] >= 4)]["delta"].mean()
                    - s[(s["ss"] == "H") & (s["idx_end"] <= 1)]["delta"].mean())
    rep["helix_middle_minus_end"] = {"difference_kcal": round(float(mv(df)), 3),
                                     "ci95_by_domain": boot_by_domain(df, mv)}

    print(json.dumps(rep, indent=2))
    json.dump(rep, open(out_path, "w"), indent=2)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
