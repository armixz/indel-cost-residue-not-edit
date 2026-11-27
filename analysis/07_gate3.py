"""Gate 3: does pLDDT beat ProDomino for short insertions on proteins that did not generate the
claim? See ../RESULTS.md. Held out: TEM-1 (enzyme fitness) and human DHFR (folding sensor).

Usage: python -I 07_gate3.py <mavedb dir> <prodomino_gate3.csv> <out.json>
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

AA1 = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Glu": "E", "Gln": "Q",
       "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F",
       "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V"}
PDB3 = {v.upper(): k for k, v in [(a, b) for a, b in zip(
    ["ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", "ILE", "LEU", "LYS", "MET",
     "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"],
    ["A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"])]}
PDB3 = dict(zip(["ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", "ILE", "LEU", "LYS",
                 "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"],
                list("ARNDCEQGHILKMFPSTWYV")))
MAXASA = dict(zip("ARNDCEQGHILKMFPSTWYV",
                  [129, 274, 195, 193, 167, 223, 225, 104, 224, 197, 201, 236, 224, 240, 159, 155, 172, 285, 263, 174]))
VDW = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80}
ORIENT = {"plddt": -1, "prodomino": +1, "rsa": +1, "is_loop": +1, "dist_term": -1}


def parse_pdb(path):
    res = {}
    for line in open(path):
        if not line.startswith("ATOM") or line[16] not in (" ", "A"):
            continue
        aa = PDB3.get(line[17:20].strip())
        if aa is None:
            continue
        i = int(line[22:26]); nm = line[12:16].strip()
        el = (line[76:78].strip() or nm[0]).upper()
        res.setdefault(i, {"aa": aa, "atoms": [], "b": []})
        res[i]["atoms"].append((nm, el, (float(line[30:38]), float(line[38:46]), float(line[46:54]))))
        res[i]["b"].append(float(line[60:66]))
    return res


def sphere(n=92):
    off, inc = 2.0 / n, math.pi * (3.0 - math.sqrt(5.0))
    pts = []
    for k in range(n):
        y = k * off - 1 + off / 2
        r = math.sqrt(max(0.0, 1 - y * y))
        pts.append([math.cos(k * inc) * r, y, math.sin(k * inc) * r])
    return np.array(pts)


def rsa_per_residue(res, probe=1.4):
    xyz, rad, own = [], [], []
    for i, d in sorted(res.items()):
        for nm, el, c in d["atoms"]:
            xyz.append(c); rad.append(VDW.get(el, 1.7) + probe); own.append(i)
    xyz, rad, own = np.array(xyz), np.array(rad), np.array(own)
    pts = sphere(); tot = np.zeros(len(xyz))
    for k in range(len(xyz)):
        d = np.linalg.norm(xyz - xyz[k], axis=1)
        nb = np.flatnonzero((d < rad + rad[k]) & (d > 0))
        test = xyz[k] + rad[k] * pts
        acc = np.ones(len(test), bool) if not len(nb) else np.all(
            np.linalg.norm(test[:, None, :] - xyz[nb][None, :, :], axis=2) >= rad[nb][None, :], axis=1)
        tot[k] = 4 * math.pi * rad[k] ** 2 * acc.mean()
    out = {}
    for i in np.unique(own):
        a = res[int(i)]["aa"]
        out[int(i)] = float(tot[own == i].sum()) / MAXASA.get(a, 200)
    return out


def dihedral(p0, p1, p2, p3):
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1n = b1 / (np.linalg.norm(b1) + 1e-9)
    v, w = b0 - np.dot(b0, b1n) * b1n, b2 - np.dot(b2, b1n) * b1n
    return math.degrees(math.atan2(np.dot(np.cross(b1n, v), w), np.dot(v, w)))


def ss_loop(res):
    get = lambda i, nm: next((np.array(a[2]) for a in res[i]["atoms"] if a[0] == nm), None) if i in res else None
    out = {}
    for i in sorted(res):
        C0, N1, CA1, C1, N2 = get(i - 1, "C"), get(i, "N"), get(i, "CA"), get(i, "C"), get(i + 1, "N")
        if any(v is None for v in (C0, N1, CA1, C1, N2)):
            out[i] = 1.0; continue
        phi, psi = dihedral(C0, N1, CA1, C1), dihedral(N1, CA1, C1, N2)
        h = -160 <= phi <= -20 and -120 <= psi <= 50
        e = -180 <= phi <= -40 and (90 <= psi <= 180 or -180 <= psi <= -150)
        out[i] = 0.0 if (h or e) else 1.0
    return out


def insertions(path):
    rows = []
    for r in pd.read_csv(path).to_dict("records"):
        h = str(r.get("hgvs_pro") or "")
        m = re.match(r"p\.([A-Za-z]{3})(\d+)_[A-Za-z]{3}\d+ins([A-Za-z]+)$", h)
        s = r.get("score")
        if not m or s is None or not np.isfinite(s):
            continue
        ins = re.findall(r"[A-Z][a-z]{2}", m.group(3))
        if len(ins) != 1:
            continue
        rows.append({"pos": int(m.group(2)), "ins_aa": AA1.get(ins[0], "X"), "score": float(s)})
    return pd.DataFrame(rows)


def cutoff(x):
    g = GaussianMixture(2, random_state=0, n_init=5).fit(np.asarray(x).reshape(-1, 1))
    mu = sorted(g.means_.ravel())
    return float(np.mean(mu)), [round(float(m), 3) for m in mu]


def evaluate(per_pos, feats, name, rep):
    cut, modes = cutoff(per_pos["score"].values)
    lab = (per_pos["score"] >= cut).astype(int)
    rep[name] = {"positions": int(len(per_pos)), "tolerated": int(lab.sum()),
                 "cutoff": round(cut, 4), "mixture_means": modes, "features": {}}
    for f, sign in ORIENT.items():
        v = per_pos[f].values.astype(float) * sign
        m = np.isfinite(v)
        if m.sum() < 20 or len(np.unique(lab.values[m])) < 2:
            continue
        rep[name]["features"][f] = {
            "auroc": round(float(roc_auc_score(lab.values[m], v[m])), 3),
            "spearman": round(float(spearmanr(v[m], per_pos["score"].values[m]).statistic), 3)}
    return rep[name]


def main(mdir, pd_path, out_path):
    pdsc = pd.read_csv(pd_path)
    rep = {}

    sets = {
        "TEM1": {"pdb": f"{mdir}/af_P62593.pdb", "pdname": "TEM1",
                 "files": [f"{mdir}/scores_00000100-a-1.csv"]},
        "DHFR": {"pdb": f"{mdir}/af_P00374.pdb", "pdname": "DHFR",
                 "files": [f"{mdir}/scores_00001219-a-{i}.csv" for i in (1, 2, 3, 4)]},
    }
    for name, cfg in sets.items():
        res = parse_pdb(cfg["pdb"])
        rsa, loop = rsa_per_residue(res), ss_loop(res)
        plddt = {i: float(np.mean(d["b"])) for i, d in res.items()}
        L = max(res)
        pdv = pdsc[pdsc["name"] == cfg["pdname"]].set_index("pos")["prodomino"]

        per_cond = []
        for f in cfg["files"]:
            d = insertions(f)
            d["cond"] = f.rsplit("-", 1)[-1].split(".")[0]
            per_cond.append(d)
        allv = pd.concat(per_cond, ignore_index=True)
        rep[name + "_raw"] = {"variants": int(len(allv)), "positions": int(allv["pos"].nunique()),
                              "conditions": int(len(per_cond))}

        # primary: TEM-1 median over inserted residues; DHFR mean over conditions
        if name == "TEM1":
            pp = allv.groupby("pos")["score"].median().reset_index()
        else:
            bycond = allv.groupby(["cond", "pos"])["score"].median().reset_index()
            pp = bycond.groupby("pos")["score"].mean().reset_index()
        pp["plddt"] = [plddt.get(p, np.nan) for p in pp["pos"]]
        pp["rsa"] = [rsa.get(p, np.nan) for p in pp["pos"]]
        pp["is_loop"] = [loop.get(p, np.nan) for p in pp["pos"]]
        pp["dist_term"] = [min(p - 1, L - p) for p in pp["pos"]]
        pp["prodomino"] = [pdv.get(p, np.nan) for p in pp["pos"]]
        evaluate(pp, None, name, rep)
        pp.to_csv(out_path.replace(".json", f"_{name}.csv"), index=False)

        if name == "TEM1":
            mature = pp[pp["pos"] > 23]
            evaluate(mature, None, "TEM1_mature_only", rep)
            # does the identity of the inserted residue matter?
            w = allv.pivot_table(index="pos", columns="ins_aa", values="score")
            cols = [c for c in w.columns if w[c].notna().sum() > 100]
            cc = w[cols].corr(method="spearman")
            iu = np.triu_indices(len(cols), 1)
            rep["TEM1_insert_identity"] = {
                "residues_compared": cols,
                "median_pairwise_spearman_between_inserted_residues": round(float(np.median(cc.values[iu])), 3),
                "range": [round(float(cc.values[iu].min()), 3), round(float(cc.values[iu].max()), 3)],
                "mean_score_by_inserted_residue": {c: round(float(w[c].mean()), 3) for c in cols}}
        else:
            for cond, g in allv.groupby("cond"):
                q = g.groupby("pos")["score"].median().reset_index()
                for c in ("plddt", "rsa", "is_loop", "dist_term", "prodomino"):
                    q[c] = pp.set_index("pos")[c].reindex(q["pos"]).values
                evaluate(q, None, f"DHFR_cond_{cond}", rep)

    verdict = {}
    for p in ("TEM1", "DHFR"):
        f = rep[p]["features"]
        verdict[p] = {"plddt": f.get("plddt", {}).get("auroc"),
                      "prodomino": f.get("prodomino", {}).get("auroc")}
        verdict[p]["plddt_at_least_prodomino"] = bool(
            verdict[p]["plddt"] is not None and verdict[p]["prodomino"] is not None
            and verdict[p]["plddt"] >= verdict[p]["prodomino"])
    verdict["supported_on_both"] = bool(all(verdict[p]["plddt_at_least_prodomino"] for p in ("TEM1", "DHFR")))
    rep["verdict"] = verdict
    json.dump(rep, open(out_path, "w"), indent=1)
    print(json.dumps(rep, indent=1))


if __name__ == "__main__":
    main(*sys.argv[1:4])
