"""Gate 10: replicate Gate 9 on the Tsuboyama domains. Marks are in ../RESULTS.md and fixed.

Usage: python -I 15_gate10_replicate.py <results dir> <out.json>
"""
import json
import sys
from collections import Counter

import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr
from sklearn.metrics import roc_auc_score

ORI = {"plddt": -1, "rsa": +1, "dist_term": -1, "prodomino": +1}
CAND = ["plddt", "rsa", "dist_term", "prodomino", "consensus"]
KS = [8, 12, 16, 24, 32]
REPEATS = 50
MIN_POS, CI_MAX, DDG_CUT = 20, 1.0, -1.0


def auroc(y, s):
    m = np.isfinite(s)
    if m.sum() < 5 or len(np.unique(y[m])) < 2:
        return np.nan
    return float(roc_auc_score(y[m], s[m]))


def main(rdir, out_path):
    ind = pd.read_csv(f"{rdir}/tsuboyama_indels.csv")
    feat = pd.read_csv(f"{rdir}/tsub_features.csv")
    wt = json.load(open(f"{rdir}/tsub_wt_dG.json"))
    pdm = pd.read_csv(f"{rdir}/prodomino_tsub.csv")

    ins = ind[(ind["kind"] == "ins") & (ind["dG_ci"] <= CI_MAX)].copy()
    ins["wt"] = ins["domain"].map(wt)
    ins = ins[ins["wt"].notna()]
    ins["ddg"] = ins["dG"] - ins["wt"]

    f = feat.set_index(["domain", "pos"])
    p = pdm.rename(columns={"name": "domain", "pos": "struct_pos"}).set_index(["domain", "struct_pos"])["prodomino"]

    rep = {"config": {"min_positions": MIN_POS, "ci_max": CI_MAX, "ddg_cut": DDG_CUT,
                      "repeats": REPEATS}, "per_k": {}}

    # ---- R4: glycine against alanine at the same positions ----
    piv = ins.pivot_table(index=["domain", "pos"], columns="res", values="ddg")
    r4 = []
    for d, g in piv.groupby(level=0):
        gg = g.dropna()
        if len(gg) >= 15:
            r4.append(float(spearmanr(gg["G"], gg["A"]).statistic))
    rep["R4"] = {"domains": len(r4), "median_spearman_G_vs_A": round(float(np.median(r4)), 3),
                 "q25": round(float(np.quantile(r4, .25)), 3), "q75": round(float(np.quantile(r4, .75)), 3),
                 "frac_above_0.8": round(float(np.mean(np.array(r4) > 0.8)), 3),
                 "passed": bool(np.median(r4) > 0.8)}
    rep["mean_ddg_by_inserted_residue"] = {r: round(float(ins[ins["res"] == r]["ddg"].mean()), 3) for r in ("G", "A")}

    # ---- build per-domain tables, pooling the two inserted residues by their mean ddG ----
    pooled = ins.groupby(["domain", "pos"])["ddg"].mean().reset_index()
    tables = {}
    for d, g in pooled.groupby("domain"):
        g = g.copy()
        try:
            ff = f.loc[d]
        except KeyError:
            continue
        g = g.join(ff, on="pos")
        g["prodomino"] = [p.get((d, sp), np.nan) for sp in g["struct_pos"]]
        g["y"] = (g["ddg"] > DDG_CUT).astype(int)
        g = g.dropna(subset=["plddt", "rsa", "dist_term", "prodomino"])
        if len(g) >= MIN_POS and 0 < g["y"].sum() < len(g):
            cols = {}
            for c in CAND:
                if c == "consensus":
                    cols[c] = np.mean([rankdata(ORI[b] * g[b].values.astype(float))
                                       for b in ("plddt", "rsa", "dist_term", "prodomino")], 0)
                else:
                    cols[c] = ORI[c] * g[c].values.astype(float)
            tables[d] = {"y": g["y"].values.astype(int), "cand": cols}
    rep["domains_used"] = len(tables)
    rep["positions_used"] = int(sum(len(t["y"]) for t in tables.values()))

    # ---- R3: which candidate is best per domain, with full knowledge ----
    best_full = Counter()
    full_auroc = {c: [] for c in CAND}
    for d, t in tables.items():
        a = {c: auroc(t["y"], t["cand"][c]) for c in CAND}
        for c in CAND:
            full_auroc[c].append(a[c])
        fin = [(v, c) for c, v in a.items() if np.isfinite(v)]
        if fin:
            best_full[max(fin)[1]] += 1
    tot = sum(best_full.values())
    rep["R3"] = {"best_candidate_counts": dict(best_full.most_common()),
                 "top_share": round(max(best_full.values()) / tot, 3) if tot else None,
                 "passed": bool(tot and max(best_full.values()) / tot <= 0.60)}
    rep["full_knowledge_auroc"] = {c: {"mean": round(float(np.nanmean(v)), 3),
                                       "p10": round(float(np.nanquantile(v, .10)), 3)}
                                   for c, v in full_auroc.items()}

    # ---- R1, R2: selection ----
    for k in KS:
        sel, fixed = [], {c: [] for c in CAND}
        oracle = []
        for d, t in tables.items():
            y, cand = t["y"], t["cand"]
            n = len(y)
            if k >= n - 8:
                continue
            s_vals, f_vals, o_vals = [], {c: [] for c in CAND}, []
            for r in range(REPEATS):
                rng = np.random.default_rng(abs(hash((d, k, r))) % (2 ** 32))
                tr = rng.choice(n, k, replace=False)
                te = np.setdiff1d(np.arange(n), tr)
                if len(np.unique(y[tr])) < 2 or len(np.unique(y[te])) < 2:
                    continue
                panel = {c: auroc(y[tr], cand[c][tr]) for c in CAND}
                fin = [(v, c) for c, v in panel.items() if np.isfinite(v)]
                if not fin:
                    continue
                s_vals.append(auroc(y[te], cand[max(fin)[1]][te]))
                for c in CAND:
                    f_vals[c].append(auroc(y[te], cand[c][te]))
                o_vals.append(max(auroc(y[te], cand[c][te]) for c in CAND))
            if s_vals:
                sel.append(np.nanmedian(s_vals))
                for c in CAND:
                    fixed[c].append(np.nanmedian(f_vals[c]))
                oracle.append(np.nanmedian(o_vals))
        if not sel:
            continue
        def stat(v):
            return {"mean": round(float(np.nanmean(v)), 3), "p10": round(float(np.nanquantile(v, .10)), 3),
                    "median": round(float(np.nanmedian(v)), 3)}
        rep["per_k"][str(k)] = {"domains": len(sel), "select": stat(sel), "oracle": stat(oracle),
                                **{f"fixed_{c}": stat(fixed[c]) for c in CAND}}
    h = rep["per_k"].get("24") or rep["per_k"][sorted(rep["per_k"])[-1]]
    rep["R1"] = {"select_mean": h["select"]["mean"], "plddt_mean": h["fixed_plddt"]["mean"],
                 "prodomino_mean": h["fixed_prodomino"]["mean"],
                 "passed": bool(h["select"]["mean"] > h["fixed_plddt"]["mean"]
                                and h["select"]["mean"] > h["fixed_prodomino"]["mean"])}
    rep["R2"] = {"select_p10": h["select"]["p10"], "plddt_p10": h["fixed_plddt"]["p10"],
                 "passed": bool(h["select"]["p10"] > h["fixed_plddt"]["p10"])}
    json.dump(rep, open(out_path, "w"), indent=1)
    print(json.dumps({k: v for k, v in rep.items() if k != "per_k"}, indent=1))
    print("\nper k:", json.dumps(rep["per_k"], indent=1))


if __name__ == "__main__":
    main(*sys.argv[1:3])
