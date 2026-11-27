"""Gate 8 and Gate 9: ESM-2 as a zero-measurement predictor, and using a small panel to CHOOSE a
predictor rather than to fit one. See ../RESULTS.md.

Usage: python -I 11_gate9_select.py <insertion-lab root> <out.json>
"""
import importlib.util
import json
import sys

import numpy as np
import pandas as pd
from scipy.stats import rankdata
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

KS = [8, 16, 24, 32, 48, 64, 96, 128]
REPEATS = 50
ORI = {"plddt": -1, "rsa": +1, "is_loop": +1, "dist_term": -1, "prodomino": +1,
       "esm_entropy": +1, "esm_wt_loglik": -1}
FIT_FEATS = ["plddt", "rsa", "is_loop", "dist_term", "prodomino"]


def auroc(y, s):
    m = np.isfinite(s)
    if m.sum() < 5 or len(np.unique(y[m])) < 2:
        return np.nan
    return float(roc_auc_score(y[m], s[m]))


def main(root, out_path):
    spec = importlib.util.spec_from_file_location("g6", f"{root}/analysis/10_gate6_how_many.py")
    g6 = importlib.util.module_from_spec(spec); spec.loader.exec_module(g6)
    data = g6.build_datasets(root)
    rep = {"gate8": {}, "gate9": {}, "selection_counts": {}}

    # ---- Gate 8: ESM as a zero-measurement predictor ----
    have_esm = all("esm_entropy" in d.columns for d in data.values())
    rep["gate8"]["available"] = bool(have_esm)
    if have_esm:
        tbl = {}
        for nm, d in data.items():
            tbl[nm] = {f: round(a, 3) for f in ("plddt", "esm_entropy", "esm_wt_loglik")
                       if np.isfinite(a := auroc(d["y"].values, ORI[f] * d[f].values.astype(float)))}
        rep["gate8"]["per_dataset"] = tbl
        for f in ("plddt", "esm_entropy", "esm_wt_loglik"):
            v = [tbl[n][f] for n in tbl if f in tbl[n]]
            rep["gate8"][f] = {"mean": round(float(np.mean(v)), 3), "worst": round(float(min(v)), 3)}
        beats = {f: sum(1 for n in tbl if tbl[n].get(f, 0) > tbl[n]["plddt"])
                 for f in ("esm_entropy", "esm_wt_loglik")}
        rep["gate8"]["beats_plddt_count"] = beats
        rep["gate8"]["passed"] = bool(max(beats.values()) >= 5)
        # decision on including ESM as a Gate 9 candidate, taken on Gate 8 numbers alone
        rep["gate8"]["included_in_gate9"] = bool(rep["gate8"]["passed"])

    CAND = ["plddt", "rsa", "dist_term", "prodomino", "consensus"]
    if have_esm and rep["gate8"].get("included_in_gate9"):
        CAND += ["esm_entropy", "esm_wt_loglik"]
    rep["gate9"]["candidates"] = CAND

    # ---- build oriented candidate matrices ----
    prepared = {}
    for nm, d in data.items():
        cols = {}
        for f in CAND:
            if f == "consensus":
                base = [c for c in ("plddt", "rsa", "dist_term", "prodomino") if c in d]
                R = np.mean([rankdata(ORI[c] * np.nan_to_num(d[c].values.astype(float),
                                                             nan=float(np.nanmedian(d[c])))) for c in base], 0)
                cols[f] = R
            else:
                cols[f] = ORI[f] * np.nan_to_num(d[f].values.astype(float), nan=float(np.nanmedian(d[f])))
        X = d[[c for c in FIT_FEATS if c in d]].astype(float)
        prepared[nm] = {"y": d["y"].values.astype(int), "cand": cols,
                        "Xfit": X.fillna(X.median()).values}

    # ---- Gate 9 ----
    results = {nm: {k: {"select": [], "fit": [], **{f"fixed_{c}": [] for c in CAND}, "oracle": []}
                    for k in KS} for nm in data}
    for nm, p in prepared.items():
        y, cand, Xf = p["y"], p["cand"], p["Xfit"]
        n = len(y)
        for k in KS:
            if k >= n - 20:
                continue
            for r in range(REPEATS):
                rng = np.random.default_rng(7000 + 100 * k + r)
                tr = rng.choice(n, k, replace=False)
                te = np.setdiff1d(np.arange(n), tr)
                if len(np.unique(y[tr])) < 2 or len(np.unique(y[te])) < 2:
                    continue
                panel = {c: auroc(y[tr], cand[c][tr]) for c in CAND}
                best = max((v, c) for c, v in panel.items() if np.isfinite(v))[1]
                results[nm][k]["select"].append(auroc(y[te], cand[best][te]))
                results[nm][k].setdefault("picked", []).append(best)
                for c in CAND:
                    results[nm][k][f"fixed_{c}"].append(auroc(y[te], cand[c][te]))
                results[nm][k]["oracle"].append(max(auroc(y[te], cand[c][te]) for c in CAND))
                m = make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000, C=1.0))
                m.fit(Xf[tr], y[tr])
                results[nm][k]["fit"].append(auroc(y[te], m.predict_proba(Xf[te])[:, 1]))

    summ = {}
    for nm in results:
        summ[nm] = {}
        for k in results[nm]:
            if not results[nm][k]["select"]:
                continue
            summ[nm][str(k)] = {m: round(float(np.nanmedian(v)), 3)
                                for m, v in results[nm][k].items() if m != "picked" and v}
        if 24 in results[nm] and results[nm][24].get("picked"):
            vals, cnt = results[nm][24]["picked"], {}
            for c in vals:
                cnt[c] = cnt.get(c, 0) + 1
            rep["selection_counts"][nm] = dict(sorted(cnt.items(), key=lambda t: -t[1]))
    rep["gate9"]["per_dataset"] = summ

    def across(method, k="24"):
        v = [summ[nm][k][method] for nm in summ if k in summ[nm] and method in summ[nm][k]]
        return (round(float(np.mean(v)), 3), round(float(min(v)), 3), len(v)) if v else (None, None, 0)

    tbl = {}
    for m in ["select", "fit", "oracle"] + [f"fixed_{c}" for c in CAND]:
        mean, worst, nds = across(m)
        tbl[m] = {"mean": mean, "worst": worst, "datasets": nds}
    rep["gate9"]["at_k24"] = tbl
    s, p = tbl["select"], tbl["fixed_plddt"]
    pdm = tbl.get("fixed_prodomino", {})
    rep["gate9"]["M3"] = {"select_mean": s["mean"], "plddt_mean": p["mean"],
                          "select_worst": s["worst"], "plddt_worst": p["worst"],
                          "passed": bool(s["mean"] and p["mean"] and s["mean"] > p["mean"]
                                         and s["worst"] > p["worst"])}
    rep["gate9"]["M4"] = {"select_worst": s["worst"], "prodomino_worst": pdm.get("worst"),
                          "passed": bool(s["worst"] and pdm.get("worst") and s["worst"] >= pdm["worst"])}
    first = None
    for k in [str(x) for x in KS]:
        a, _, _ = across("select", k); b, _, _ = across("fixed_plddt", k)
        if a and b and a > b:
            first = int(k); break
    rep["gate9"]["k_where_selection_overtakes_plddt"] = first
    json.dump(rep, open(out_path, "w"), indent=1)
    print(json.dumps({k: v for k, v in rep.items() if k != "gate9"}, indent=1))
    print("\ngate9 at k=24:", json.dumps(tbl, indent=1))
    print("M3:", rep["gate9"]["M3"]); print("M4:", rep["gate9"]["M4"])
    print("overtakes pLDDT at k =", first)
    print("\nper dataset, select vs fixed_plddt vs oracle:")
    for nm in summ:
        if "24" in summ[nm]:
            r = summ[nm]["24"]
            print(f"  {nm:16s} select={r.get('select')} plddt={r.get('fixed_plddt')} "
                  f"prodomino={r.get('fixed_prodomino')} oracle={r.get('oracle')} fit={r.get('fit')}")


if __name__ == "__main__":
    main(*sys.argv[1:3])
