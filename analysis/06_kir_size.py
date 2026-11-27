"""Gate 2: does insertion tolerance depend on insert size, and does a domain-insertion model
transfer to short inserts? See ../RESULTS.md.

Kir2.1 (DIMPLE, MaveDB urn:mavedb:00000660-a-1): surface-expression scores for insertions of
Gly (1 aa), GlySer (2 aa) and GlySerGly (3 aa) at the same positions.

Usage: python -I 06_kir_size.py <kir21_scores.csv> <af_P35561.pdb> <prodomino_kir21.csv> <out.json>
"""
import json
import math
import re
import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score
from sklearn.mixture import GaussianMixture

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from importlib import import_module
F = import_module("01_features") if False else None   # structure helpers are re-implemented below

AA3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLU": "E", "GLN": "Q",
       "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
       "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}
MAXASA = {"A": 129, "R": 274, "N": 195, "D": 193, "C": 167, "E": 223, "Q": 225, "G": 104, "H": 224,
          "I": 197, "L": 201, "K": 236, "M": 224, "F": 240, "P": 159, "S": 155, "T": 172, "W": 285,
          "Y": 263, "V": 174}
VDW = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80}
FLAG_AFTER, FLAG_LEN = 115, 8          # DYKDDDDK inserted after construct position 115


def to_struct(p):
    """construct numbering -> AlphaFold residue number, or None for the tag itself."""
    if p <= FLAG_AFTER:
        return p
    if p <= FLAG_AFTER + FLAG_LEN:
        return None
    return p - FLAG_LEN


def parse_pdb(path):
    res = {}
    for line in open(path):
        if not line.startswith("ATOM") or line[16] not in (" ", "A"):
            continue
        rn = AA3.get(line[17:20].strip())
        if rn is None:
            continue
        i = int(line[22:26])
        nm = line[12:16].strip()
        el = (line[76:78].strip() or nm[0]).upper()
        res.setdefault(i, {"aa": rn, "atoms": [], "b": []})
        res[i]["atoms"].append((nm, el, (float(line[30:38]), float(line[38:46]), float(line[46:54]))))
        res[i]["b"].append(float(line[60:66]))
    return res


def sphere(n=92):
    off, inc = 2.0 / n, math.pi * (3.0 - math.sqrt(5.0))
    return np.array([[math.cos(k * inc) * math.sqrt(max(0.0, 1 - (k * off - 1 + off / 2) ** 2)),
                      k * off - 1 + off / 2,
                      math.sin(k * inc) * math.sqrt(max(0.0, 1 - (k * off - 1 + off / 2) ** 2))]
                     for k in range(n)])


def sasa_per_residue(res, probe=1.4):
    xyz, rad, own = [], [], []
    for i, d in sorted(res.items()):
        for nm, el, c in d["atoms"]:
            xyz.append(c); rad.append(VDW.get(el, 1.7) + probe); own.append(i)
    xyz, rad, own = np.array(xyz), np.array(rad), np.array(own)
    pts, out = sphere(), {}
    tot = np.zeros(len(xyz))
    for k in range(len(xyz)):
        d = np.linalg.norm(xyz - xyz[k], axis=1)
        nb = np.flatnonzero((d < rad + rad[k]) & (d > 0))
        test = xyz[k] + rad[k] * pts
        acc = np.ones(len(test), bool) if not len(nb) else \
            np.all(np.linalg.norm(test[:, None, :] - xyz[nb][None, :, :], axis=2) >= rad[nb][None, :], axis=1)
        tot[k] = 4 * math.pi * rad[k] ** 2 * acc.mean()
    for i in np.unique(own):
        out[int(i)] = float(tot[own == i].sum())
    return out


def dihedral(p0, p1, p2, p3):
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1n = b1 / (np.linalg.norm(b1) + 1e-9)
    v, w = b0 - np.dot(b0, b1n) * b1n, b2 - np.dot(b2, b1n) * b1n
    return math.degrees(math.atan2(np.dot(np.cross(b1n, v), w), np.dot(v, w)))


def ss_from_dihedrals(res):
    get = lambda i, nm: next((np.array(a[2]) for a in res[i]["atoms"] if a[0] == nm), None) if i in res else None
    out = {}
    for i in sorted(res):
        C0, N1, CA1, C1, N2 = get(i - 1, "C"), get(i, "N"), get(i, "CA"), get(i, "C"), get(i + 1, "N")
        if any(v is None for v in (C0, N1, CA1, C1, N2)):
            out[i] = "L"; continue
        phi, psi = dihedral(C0, N1, CA1, C1), dihedral(N1, CA1, C1, N2)
        out[i] = "H" if (-160 <= phi <= -20 and -120 <= psi <= 50) else \
                 ("E" if (-180 <= phi <= -40 and (90 <= psi <= 180 or -180 <= psi <= -150)) else "L")
    return out


def load_scores(path):
    df = pd.read_csv(path)
    rows = []
    for p, s in zip(df["hgvs_pro"].fillna(""), df["score"]):
        if not np.isfinite(s):
            continue
        m = re.match(r"p\.([A-Za-z]{3})(\d+)_([A-Za-z]{3})(\d+)ins([A-Za-z]+)$", p)
        if m:
            n = len(re.findall(r"[A-Z][a-z]{2}", m.group(5)))
            rows.append({"pos": int(m.group(2)), "kind": f"ins{n}", "score": float(s)})
            continue
        d = re.match(r"p\.([A-Za-z]{3})(\d+)del$", p)
        if d:
            rows.append({"pos": int(d.group(2)), "kind": "del1", "score": float(s)})
    return pd.DataFrame(rows)


def mixture_cutoff(x):
    g = GaussianMixture(2, random_state=0, n_init=5).fit(np.asarray(x).reshape(-1, 1))
    mu = sorted(g.means_.ravel())
    return float(np.mean(mu)), [round(float(m), 3) for m in mu]


def boot_spear(a, b, n=2000, seed=0):
    rng = np.random.default_rng(seed)
    v = [spearmanr(a[i], b[i]).statistic for i in (rng.integers(0, len(a), len(a)) for _ in range(n))]
    return [round(float(q), 3) for q in np.nanquantile(v, [0.025, 0.975])]


def main(scores_path, pdb_path, pd_path, out_path):
    sc = load_scores(scores_path)
    wide = sc.pivot_table(index="pos", columns="kind", values="score")
    rep = {"n_positions": int(len(wide)),
           "coverage": {k: int(wide[k].notna().sum()) for k in wide.columns}}

    # Q1: size dependence
    rep["Q1_size_dependence"] = {}
    for a, b in (("ins1", "ins2"), ("ins2", "ins3"), ("ins1", "ins3"), ("ins1", "del1"), ("ins3", "del1")):
        m = wide[a].notna() & wide[b].notna()
        x, y = wide.loc[m, a].values, wide.loc[m, b].values
        rep["Q1_size_dependence"][f"{a}_vs_{b}"] = {
            "n": int(m.sum()), "spearman": round(float(spearmanr(x, y).statistic), 3),
            "ci95": boot_spear(x, y), "pearson": round(float(np.corrcoef(x, y)[0, 1]), 3)}
    rep["median_scores"] = {k: round(float(wide[k].median()), 3) for k in wide.columns}

    # threshold from the measured scores only
    allins = np.concatenate([wide[k].dropna().values for k in ("ins1", "ins2", "ins3") if k in wide])
    cut, modes = mixture_cutoff(allins)
    rep["mixture_cutoff"] = {"cutoff": round(cut, 3), "component_means": modes,
                             "fraction_below": {k: round(float((wide[k] < cut).mean()), 3) for k in wide.columns}}

    # features
    res = parse_pdb(pdb_path)
    sasa, ss = sasa_per_residue(res), ss_from_dihedrals(res)
    plddt = {i: float(np.mean(d["b"])) for i, d in res.items()}
    L = len(wide.index.max() and wide.index) and int(wide.index.max())
    feat = {}
    for p in wide.index:
        s = to_struct(int(p))
        aa = res.get(s, {}).get("aa") if s else None
        feat[p] = {
            "rsa": (sasa.get(s) / MAXASA[aa]) if (s and aa in MAXASA and sasa.get(s) is not None) else np.nan,
            "plddt": plddt.get(s, np.nan) if s else np.nan,
            "is_loop": float(ss.get(s) == "L") if s else np.nan,
            "dist_term": min(int(p) - 1, int(wide.index.max()) - int(p)),
            "in_tag": s is None}
    fdf = pd.DataFrame(feat).T

    pdsc = pd.read_csv(pd_path).set_index("pos")["prodomino"]
    fdf["prodomino"] = [pdsc.get(p, np.nan) for p in fdf.index]

    ORIENT = {"prodomino": +1, "rsa": +1, "plddt": -1, "is_loop": +1, "dist_term": -1}
    rep["Q2_Q3"] = {}
    for kind in ("ins1", "ins2", "ins3", "del1"):
        if kind not in wide:
            continue
        y = wide[kind]
        ok = y.notna() & ~fdf["in_tag"].astype(bool)
        lab = (y[ok] >= cut).astype(int)                 # 1 = tolerated
        rep["Q2_Q3"][kind] = {"n": int(ok.sum()), "tolerated": int(lab.sum())}
        for f, sign in ORIENT.items():
            v = pd.to_numeric(fdf.loc[ok, f], errors="coerce").values * sign
            m = np.isfinite(v)
            if m.sum() < 20 or len(np.unique(lab.values[m])) < 2:
                continue
            rep["Q2_Q3"][kind][f] = {
                "spearman_vs_score": round(float(spearmanr(v[m], y[ok].values[m]).statistic), 3),
                "auroc_tolerated": round(float(roc_auc_score(lab.values[m], v[m])), 3)}
    pdm = rep["Q2_Q3"].get("ins1", {}).get("prodomino", {})
    rep["verdict_Q2"] = {"prodomino_auroc_ins1": pdm.get("auroc_tolerated"),
                         "pass_mark_0.65": bool((pdm.get("auroc_tolerated") or 0) >= 0.65)}
    fdf.join(wide).to_csv(out_path.replace(".json", "_table.csv"))
    json.dump(rep, open(out_path, "w"), indent=1)
    print(json.dumps(rep, indent=1))


if __name__ == "__main__":
    main(*sys.argv[1:5])
