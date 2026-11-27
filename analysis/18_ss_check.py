"""Independent validation of the DSSP-style assignment in 17_dssp.py.

Checks the hydrogen-bond assignment against backbone geometry it never used:
  - residues called H should centre near phi -63, psi -42
  - residues called E should centre near phi -120, psi +130
  - helical residues should show CA(i)-CA(i+4) near 6.2 A, strands near 13 A

Usage: run.sh 18_ss_check.py <afdir> [nmodels]
"""
import glob
import importlib.util
import math
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("d17", os.path.join(HERE, "17_dssp.py"))
d17 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(d17)


def dih(p0, p1, p2, p3):
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1n = b1 / (np.linalg.norm(b1) + 1e-9)
    v, w = b0 - np.dot(b0, b1n) * b1n, b2 - np.dot(b2, b1n) * b1n
    return math.degrees(math.atan2(np.dot(np.cross(b1n, v), w), np.dot(v, w)))


def main(afdir, n=250):
    paths = sorted(glob.glob(os.path.join(afdir, "*.pdb")))[:int(n)]
    rows = []
    for p in paths:
        res = d17.parse_backbone(p)
        ss, keys = d17.assign(res)
        if not ss:
            continue
        el = d17.elements(ss, keys)
        for i in keys:
            if (i - 1) not in res or (i + 1) not in res:
                continue
            phi = dih(res[i - 1]["C"], res[i]["N"], res[i]["CA"], res[i]["C"])
            psi = dih(res[i]["N"], res[i]["CA"], res[i]["C"], res[i + 1]["N"])
            d4 = (np.linalg.norm(res[i + 4]["CA"] - res[i]["CA"])
                  if (i + 4) in res else np.nan)
            rows.append((el[i][0], phi, psi, d4))
    df = pd.DataFrame(rows, columns=["ss", "phi", "psi", "ca_i_i4"])
    out = df.groupby("ss").agg(n=("phi", "size"), phi_med=("phi", "median"),
                               psi_med=("psi", "median"),
                               ca4_med=("ca_i_i4", "median")).round(1)
    print(out)
    print("\nexpected: H near phi -63 / psi -42 / CA(i,i+4) ~6.2 A;"
          " E near phi -120 / psi +130 / CA(i,i+4) ~13 A")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else 250)
