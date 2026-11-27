"""Kabsch-Sander hydrogen-bond secondary structure on the Tsuboyama AlphaFold models.

Prerequisite step of PREREG_indel.md. Produces, per domain and residue:
  ss        H / E / L by DSSP rules
  elem_id   index of the secondary-structure element the residue belongs to
  elem_len  its length
  idx_n     distance (residues) to the element's N-terminal end, 0 at the end itself
  idx_c     distance to the element's C-terminal end
  idx_end   min(idx_n, idx_c) -- distance to the NEAREST element terminus
  d_loop    distance in residues to the nearest loop (L) residue, 0 if itself loop

No ddG is read here. Usage: run.sh 17_dssp.py
"""
import glob
import json
import math
import os
import sys
from collections import defaultdict

import numpy as np

PDB3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
        "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
        "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}

Q1Q2F = 0.084 * 332.0  # kcal/mol, Kabsch & Sander 1983
HB_CUT = -0.5


def parse_backbone(path):
    """Return {resnum: {'aa':..., 'N':xyz, 'CA':xyz, 'C':xyz, 'O':xyz}} for one chain."""
    res = {}
    for line in open(path):
        if not line.startswith("ATOM") or line[16] not in (" ", "A"):
            continue
        aa = PDB3.get(line[17:20].strip())
        if aa is None:
            continue
        nm = line[12:16].strip()
        if nm not in ("N", "CA", "C", "O"):
            continue
        i = int(line[22:26])
        res.setdefault(i, {"aa": aa})
        res[i][nm] = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return {i: r for i, r in res.items() if all(k in r for k in ("N", "CA", "C", "O"))}


def amide_h(res, keys):
    """DSSP rule: H sits 1.0 A from N, opposite the preceding residue's C=O direction."""
    h = {}
    for k, i in enumerate(keys):
        if k == 0:
            continue
        p = keys[k - 1]
        if p != i - 1:
            continue  # chain break: no amide hydrogen defined
        if res[i]["aa"] == "P":
            continue  # proline has no amide hydrogen
        co = res[p]["C"] - res[p]["O"]
        n = np.linalg.norm(co)
        if n < 1e-6:
            continue
        h[i] = res[i]["N"] + co / n
    return h


def hbonds(res, keys):
    """Set of (donor, acceptor) residue pairs: donor N-H ... O=C acceptor."""
    h = amide_h(res, keys)
    ca = np.array([res[i]["CA"] for i in keys])
    out = set()
    energies = {}
    for a, i in enumerate(keys):
        if i not in h:
            continue
        d = ca - ca[a]
        near = np.where((d * d).sum(1) < 81.0)[0]  # CA within 9 A
        for b in near:
            j = keys[b]
            if abs(j - i) < 2:
                continue
            rON = np.linalg.norm(res[j]["O"] - res[i]["N"])
            rCH = np.linalg.norm(res[j]["C"] - h[i])
            rOH = np.linalg.norm(res[j]["O"] - h[i])
            rCN = np.linalg.norm(res[j]["C"] - res[i]["N"])
            if min(rON, rCH, rOH, rCN) < 0.5:
                continue
            e = Q1Q2F * (1.0 / rON + 1.0 / rCH - 1.0 / rOH - 1.0 / rCN)
            if e < HB_CUT:
                out.add((i, j))
                energies[(i, j)] = e
    return out, energies


def assign(res):
    keys = sorted(res)
    if len(keys) < 5:
        return {}, keys
    hb, _ = hbonds(res, keys)
    HB = lambda i, j: (i, j) in hb

    ss = {i: "L" for i in keys}

    # alpha helix: two consecutive 4-turns. 4-turn at i means HB(i+4 -> i).
    turn4 = {i: HB(i + 4, i) for i in keys}
    for i in keys:
        if turn4.get(i) and turn4.get(i - 1):
            for k in range(i, i + 4):
                if k in ss:
                    ss[k] = "H"

    # beta bridges
    bridge = defaultdict(set)
    for a, i in enumerate(keys):
        for j in keys:
            if j <= i + 2:
                continue
            par = (HB(i, j - 1) and HB(j + 1, i)) or (HB(j, i - 1) and HB(i + 1, j))
            anti = (HB(i, j) and HB(j, i)) or (HB(i - 1, j + 1) and HB(j - 1, i + 1))
            if par or anti:
                bridge[i].add(j)
                bridge[j].add(i)
    # a residue is E if it is in a ladder, i.e. a bridge that continues to a neighbour
    for i in keys:
        if not bridge[i]:
            continue
        cont = any((j + 1 in bridge.get(i + 1, ())) or (j - 1 in bridge.get(i + 1, ())) or
                   (j + 1 in bridge.get(i - 1, ())) or (j - 1 in bridge.get(i - 1, ()))
                   for j in bridge[i])
        if cont and ss[i] != "H":
            ss[i] = "E"
    return ss, keys


def elements(ss, keys, min_h=5, min_e=2):
    """Maximal runs of the same state. Short helices/strands are demoted to L."""
    out = []
    run = []
    for i in keys:
        if run and ss[i] == ss[run[-1]] and i == run[-1] + 1:
            run.append(i)
        else:
            if run:
                out.append(run)
            run = [i]
    if run:
        out.append(run)
    final = {}
    eid = 0
    for run in out:
        s = ss[run[0]]
        if s == "H" and len(run) < min_h:
            s = "L"
        if s == "E" and len(run) < min_e:
            s = "L"
        if s == "L":
            for i in run:
                final[i] = ("L", -1, len(run), -1, -1)
            continue
        for k, i in enumerate(run):
            final[i] = (s, eid, len(run), k, len(run) - 1 - k)
        eid += 1
    return final


def main(afdir, out_prefix):
    paths = sorted(glob.glob(os.path.join(afdir, "**", "*.pdb"), recursive=True))
    paths = [p for p in paths if "__MACOSX" not in p]
    print(f"{len(paths)} model files", flush=True)

    rows = []
    agree_num = agree_den = 0
    nfail = 0
    for n, p in enumerate(paths):
        dom = os.path.basename(p)
        try:
            res = parse_backbone(p)
            ss, keys = assign(res)
            if not ss:
                nfail += 1
                continue
            el = elements(ss, keys)
        except Exception as e:  # noqa: BLE001
            nfail += 1
            continue
        loops = [i for i in keys if el[i][0] == "L"]
        for i in keys:
            s, eid, elen, idx_n, idx_c = el[i]
            dl = 0 if s == "L" else (min(abs(i - j) for j in loops) if loops else -1)
            rows.append({"domain": dom, "struct_pos": i, "aa": res[i]["aa"], "ss": s,
                         "elem_id": eid, "elem_len": elen, "idx_n": idx_n, "idx_c": idx_c,
                         "idx_end": (min(idx_n, idx_c) if s != "L" else -1), "d_loop": dl})
        if n % 100 == 0:
            print(f"  {n}/{len(paths)}", flush=True)

    import pandas as pd
    df = pd.DataFrame(rows)
    df.to_csv(out_prefix + "_ss.csv", index=False)

    def elem_stats(code):
        sub = df[df["ss"] == code]
        if len(sub) == 0:
            return 0, None
        per = sub.groupby(["domain", "elem_id"])["elem_len"].first()
        return len(per), float(per.median())

    nh, mh = elem_stats("H")
    ne, me = elem_stats("E")
    rep = {"models": len(paths), "failed": nfail, "domains": int(df["domain"].nunique()),
           "residues": int(len(df)),
           "ss_fractions": {str(k): round(float(v), 4) for k, v in
                            df["ss"].value_counts(normalize=True).to_dict().items()},
           "helix_elements": nh, "median_helix_len": mh,
           "strand_elements": ne, "median_strand_len": me}
    print(json.dumps(rep, indent=2))
    json.dump(rep, open(out_prefix + "_ss.json", "w"), indent=2)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
