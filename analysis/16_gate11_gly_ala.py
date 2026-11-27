"""Gate 11: does the preferred inserted residue depend on local secondary structure?
delta = ddG(Ala) - ddG(Gly) at the same position. Marks are fixed in ../RESULTS.md.

Usage: python -I 16_gate11_gly_ala.py <results dir> <AF dir> <out.json>
"""
import json
import math
import sys

import numpy as np
import pandas as pd

PDB3 = dict(zip(["ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", "ILE", "LEU", "LYS",
                 "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"], list("ARNDCEQGHILKMFPSTWYV")))
CI_MAX = 1.0


def backbone(path):
    res = {}
    for line in open(path):
        if not line.startswith("ATOM") or line[16] not in (" ", "A"):
            continue
        aa = PDB3.get(line[17:20].strip())
        if aa is None:
            continue
        i, nm = int(line[22:26]), line[12:16].strip()
        if nm in ("N", "CA", "C"):
            res.setdefault(i, {"aa": aa})[nm] = np.array(
                [float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return res


def dih(p0, p1, p2, p3):
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1n = b1 / (np.linalg.norm(b1) + 1e-9)
    v, w = b0 - np.dot(b0, b1n) * b1n, b2 - np.dot(b2, b1n) * b1n
    return math.degrees(math.atan2(np.dot(np.cross(b1n, v), w), np.dot(v, w)))


def ss3(res):
    out, phi_map = {}, {}
    for i in sorted(res):
        a, b, c = res.get(i - 1), res.get(i), res.get(i + 1)
        if not (a and b and c) or any(k not in a for k in ("C",)) or any(
                k not in b for k in ("N", "CA", "C")) or "N" not in c:
            out[i] = "L"; continue
        phi = dih(a["C"], b["N"], b["CA"], b["C"])
        psi = dih(b["N"], b["CA"], b["C"], c["N"])
        phi_map[i] = phi
        if -160 <= phi <= -20 and -120 <= psi <= 50:
            out[i] = "H"
        elif -180 <= phi <= -40 and (90 <= psi <= 180 or -180 <= psi <= -150):
            out[i] = "E"
        else:
            out[i] = "L"
    return out, phi_map


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


def main(rdir, afdir, out_path):
    ind = pd.read_csv(f"{rdir}/tsuboyama_indels.csv")
    feat = pd.read_csv(f"{rdir}/tsub_features.csv")
    wt = json.load(open(f"{rdir}/tsub_wt_dG.json"))
    ins = ind[(ind["kind"] == "ins") & (ind["dG_ci"] <= CI_MAX)].copy()
    ins["wt"] = ins["domain"].map(wt)
    ins = ins[ins["wt"].notna()]
    ins["ddg"] = ins["dG"] - ins["wt"]
    piv = ins.pivot_table(index=["domain", "pos"], columns="res", values="ddg").dropna().reset_index()
    piv["delta"] = piv["A"] - piv["G"]

    # structure per domain
    ss_rows = []
    for d in piv["domain"].unique():
        try:
            res = backbone(f"{afdir}/{d}.pdb")
        except FileNotFoundError:
            continue
        s, phi = ss3(res)
        for i in sorted(res):
            ss_rows.append({"domain": d, "struct_pos": i, "ss": s.get(i, "L"),
                            "phi": phi.get(i, np.nan)})
    ss = pd.DataFrame(ss_rows)
    f = feat[["domain", "pos", "struct_pos", "rsa", "length"]]
    df = piv.merge(f, on=["domain", "pos"]).merge(ss, on=["domain", "struct_pos"])
    # flanking positive phi: glycine should be favoured most where neighbours need left-handed phi
    pp = ss.set_index(["domain", "struct_pos"])["phi"]
    df["phi_prev"] = [pp.get((d, s - 1), np.nan) for d, s in zip(df["domain"], df["struct_pos"])]
    df["phi_next"] = [pp.get((d, s + 1), np.nan) for d, s in zip(df["domain"], df["struct_pos"])]
    df["flank_pos_phi"] = ((df["phi_prev"] > 0) | (df["phi_next"] > 0)).astype(int)

    rep = {"positions": int(len(df)), "domains": int(df["domain"].nunique()),
           "mean_ddg": {"G": round(float(piv["G"].mean()), 3), "A": round(float(piv["A"].mean()), 3)},
           "overall_delta_mean": round(float(df["delta"].mean()), 3)}

    by = df.groupby("ss")["delta"].agg(["count", "mean", "median"])
    rep["by_secondary_structure"] = {k: {"n": int(v["count"]), "mean_delta": round(float(v["mean"]), 3),
                                         "median_delta": round(float(v["median"]), 3)}
                                     for k, v in by.iterrows()}
    hl = lambda s: (s[s["ss"] == "H"]["delta"].mean() - s[s["ss"] == "L"]["delta"].mean())
    diff = float(hl(df))
    ci = boot_by_domain(df, hl)
    rep["A_helix_minus_loop"] = {"difference_kcal": round(diff, 3), "ci95_by_domain": ci,
                                 "passed": bool(ci[0] > 0)}
    rep["B_effect_at_least_0.3"] = {"difference_kcal": round(diff, 3), "passed": bool(diff >= 0.3)}
    rep["H11_supported"] = bool(ci[0] > 0 and diff >= 0.3)

    # secondary analyses
    q = df["rsa"].quantile([0.33, 0.66]).values
    df["burial"] = np.where(df["rsa"] <= q[0], "buried", np.where(df["rsa"] >= q[1], "exposed", "mid"))
    rep["by_burial"] = {k: {"n": int(len(g)), "mean_delta": round(float(g["delta"].mean()), 3)}
                        for k, g in df.groupby("burial")}
    rep["by_flanking_positive_phi"] = {
        str(k): {"n": int(len(g)), "mean_delta": round(float(g["delta"].mean()), 3)}
        for k, g in df.groupby("flank_pos_phi")}
    fp = lambda s: (s[s["flank_pos_phi"] == 1]["delta"].mean() - s[s["flank_pos_phi"] == 0]["delta"].mean())
    rep["flank_phi_difference"] = {"difference_kcal": round(float(fp(df)), 3),
                                   "ci95_by_domain": boot_by_domain(df, fp)}
    # helix middle against helix ends
    df["ss_prev"] = [None] * len(df)
    ssmap = ss.set_index(["domain", "struct_pos"])["ss"]
    mid = []
    for d, s, k in zip(df["domain"], df["struct_pos"], df["ss"]):
        if k != "H":
            mid.append(np.nan); continue
        a, b = ssmap.get((d, s - 1)), ssmap.get((d, s + 1))
        mid.append(1.0 if (a == "H" and b == "H") else 0.0)
    df["helix_middle"] = mid
    h = df[df["ss"] == "H"]
    rep["helix_middle_vs_end"] = {
        str(k): {"n": int(len(g)), "mean_delta": round(float(g["delta"].mean()), 3)}
        for k, g in h.groupby("helix_middle")}
    df.to_csv(out_path.replace(".json", "_positions.csv"), index=False)
    json.dump(rep, open(out_path, "w"), indent=1)
    print(json.dumps(rep, indent=1))


if __name__ == "__main__":
    main(*sys.argv[1:4])
