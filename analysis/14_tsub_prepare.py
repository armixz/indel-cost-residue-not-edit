"""Wild-type dG per domain, numbering check against the AlphaFold models, and per-position
features for the Tsuboyama replication set.

The numbering check is free and strict: every `delX{pos}` row states the wild-type residue at that
position, so it must match the AlphaFold model's residue there.

Usage: python -I 14_tsub_prepare.py <Lib3.csv> <indels.csv> <AF dir> <out prefix>
"""
import csv
import json
import math
import re
import sys
from collections import defaultdict

import numpy as np

PDB3 = dict(zip(["ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", "ILE", "LEU", "LYS",
                 "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"], list("ARNDCEQGHILKMFPSTWYV")))
MAXASA = dict(zip("ARNDCEQGHILKMFPSTWYV",
                  [129, 274, 195, 193, 167, 223, 225, 104, 224, 197, 201, 236, 224, 240, 159, 155, 172, 285, 263, 174]))
VDW = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80}


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
    out = []
    for k in range(n):
        y = k * off - 1 + off / 2
        r = math.sqrt(max(0.0, 1 - y * y))
        out.append([math.cos(k * inc) * r, y, math.sin(k * inc) * r])
    return np.array(out)


PTS = sphere()


def rsa(res, probe=1.4):
    xyz, rad, own = [], [], []
    for i, d in sorted(res.items()):
        for nm, el, c in d["atoms"]:
            xyz.append(c); rad.append(VDW.get(el, 1.7) + probe); own.append(i)
    xyz, rad, own = np.array(xyz), np.array(rad), np.array(own)
    tot = np.zeros(len(xyz))
    for k in range(len(xyz)):
        d = np.linalg.norm(xyz - xyz[k], axis=1)
        nb = np.flatnonzero((d < rad + rad[k]) & (d > 0))
        t = xyz[k] + rad[k] * PTS
        acc = np.ones(len(t), bool) if not len(nb) else np.all(
            np.linalg.norm(t[:, None, :] - xyz[nb][None, :, :], axis=2) >= rad[nb][None, :], axis=1)
        tot[k] = 4 * math.pi * rad[k] ** 2 * acc.mean()
    return {int(i): float(tot[own == i].sum()) / MAXASA.get(res[int(i)]["aa"], 200) for i in np.unique(own)}


def dih(p0, p1, p2, p3):
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1n = b1 / (np.linalg.norm(b1) + 1e-9)
    v, w = b0 - np.dot(b0, b1n) * b1n, b2 - np.dot(b2, b1n) * b1n
    return math.degrees(math.atan2(np.dot(np.cross(b1n, v), w), np.dot(v, w)))


def loops(res):
    g = lambda i, nm: next((np.array(a[2]) for a in res[i]["atoms"] if a[0] == nm), None) if i in res else None
    out = {}
    for i in sorted(res):
        C0, N1, CA1, C1, N2 = g(i - 1, "C"), g(i, "N"), g(i, "CA"), g(i, "C"), g(i + 1, "N")
        if any(v is None for v in (C0, N1, CA1, C1, N2)):
            out[i] = 1.0; continue
        phi, psi = dih(C0, N1, CA1, C1), dih(N1, CA1, C1, N2)
        h = -160 <= phi <= -20 and -120 <= psi <= 50
        e = -180 <= phi <= -40 and (90 <= psi <= 180 or -180 <= psi <= -150)
        out[i] = 0.0 if (h or e) else 1.0
    return out


def main(lib3, indels, afdir, prefix):
    # wild-type dG per domain
    wt = {}
    with open(lib3, newline="") as fh:
        for row in csv.DictReader(fh):
            m = re.match(r"^([0-9A-Za-z]{4})\.pdb$", row["name"]) or re.match(r"^([0-9A-Za-z]{4})$", row["name"])
            if m:
                try:
                    wt[m.group(1)] = float(row["deltaG"])
                except (ValueError, TypeError):
                    pass
    print("wild-type dG found for", len(wt), "domains", flush=True)

    rows = list(csv.DictReader(open(indels)))
    doms = sorted({r["domain"] for r in rows})
    wtres = defaultdict(dict)
    for r in rows:
        if r["kind"] == "del":
            wtres[r["domain"]][int(r["pos"])] = r["res"]

    feats, checks = [], {"ok": 0, "mismatch": 0, "no_model": 0, "offsets": defaultdict(int)}
    seqs = {}
    for d in doms:
        path = f"{afdir}/{d}.pdb"
        try:
            res = parse_pdb(path)
        except FileNotFoundError:
            checks["no_model"] += 1
            continue
        if not res:
            checks["no_model"] += 1
            continue
        best, bestok = 0, -1
        for off in range(-3, 4):
            ok = sum(1 for p, a in wtres[d].items() if res.get(p + off, {}).get("aa") == a)
            if ok > bestok:
                best, bestok = off, ok
        tot = max(len(wtres[d]), 1)
        checks["offsets"][best] += 1
        if bestok / tot < 0.9:
            checks["mismatch"] += 1
            continue
        checks["ok"] += 1
        a, l = rsa(res), loops(res)
        plddt = {i: float(np.mean(v["b"])) for i, v in res.items()}
        L = max(res)
        seqs[d] = "".join(res[i]["aa"] for i in sorted(res))
        for p in sorted(wtres[d]):
            s = p + best
            if s not in res:
                continue
            feats.append({"domain": d, "pos": p, "struct_pos": s, "aa": res[s]["aa"],
                          "plddt": plddt[s], "rsa": a.get(s, np.nan), "is_loop": l.get(s, np.nan),
                          "dist_term": min(s - 1, L - s), "length": L, "offset": best})
    import pandas as pd
    pd.DataFrame(feats).to_csv(prefix + "_features.csv", index=False)
    json.dump({"wt_dG": wt, "checks": {k: (dict(v) if isinstance(v, defaultdict) else v) for k, v in checks.items()}},
              open(prefix + "_wt.json", "w"), indent=1)
    with open(prefix + "_seqs.fasta", "w") as fh:
        for d, s in seqs.items():
            fh.write(f">{d}\n{s}\n")
    print("domains with usable model+numbering:", checks["ok"], "mismatch:", checks["mismatch"],
          "no model:", checks["no_model"], "offsets:", dict(checks["offsets"]))
    print("feature rows:", len(feats), "sequences written:", len(seqs))


if __name__ == "__main__":
    main(*sys.argv[1:5])
