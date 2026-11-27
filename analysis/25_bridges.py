"""Beta-bridge partners, parallel/antiparallel type, and edge-vs-interior, per strand residue.

Needed by marks M1-M3 of PREREG_bulge.md. Uses the same Kabsch-Sander hydrogen bonds as
17_dssp.py and records what that script discards.

Per strand residue:
  n_partners   1 = edge strand position, 2 = interior
  bridge_type  'anti', 'par' or 'mixed'
  has_bulge    the wild-type structure already shows a bulge at this position: a residue whose
               bridge partner is shared with a sequence neighbour (two residues of one strand
               facing one of the other)

Usage: run.sh 25_bridges.py <afdir> <out_csv>
"""
import glob
import importlib.util
import json
import os
import sys
from collections import defaultdict

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("d17", os.path.join(HERE, "17_dssp.py"))
d17 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(d17)


def bridges(res, keys, hb):
    """Return {i: [(j, 'anti'|'par'), ...]} using the DSSP bridge definitions."""
    HB = lambda a, b: (a, b) in hb
    out = defaultdict(list)
    for i in keys:
        for j in keys:
            if j <= i + 2:
                continue
            par = (HB(i, j - 1) and HB(j + 1, i)) or (HB(j, i - 1) and HB(i + 1, j))
            anti = (HB(i, j) and HB(j, i)) or (HB(i - 1, j + 1) and HB(j - 1, i + 1))
            if anti:
                out[i].append((j, "anti"))
                out[j].append((i, "anti"))
            elif par:
                out[i].append((j, "par"))
                out[j].append((i, "par"))
    return out


def main(afdir, out_csv):
    paths = sorted(glob.glob(os.path.join(afdir, "*.pdb")))
    rows = []
    for n, p in enumerate(paths):
        dom = os.path.basename(p).replace(".pdb", "")
        res = d17.parse_backbone(p)
        ss, keys = d17.assign(res)
        if not ss:
            continue
        el = d17.elements(ss, keys)
        hb, _ = d17.hbonds(res, keys)
        br = bridges(res, keys, hb)
        for i in keys:
            if el[i][0] != "E":
                continue
            pr = br.get(i, [])
            partners = sorted({j for j, _ in pr})
            types = {t for _, t in pr}
            # a bulge: this residue and a sequence neighbour share a bridge partner
            shared = False
            for nb in (i - 1, i + 1):
                if nb in br:
                    if set(j for j, _ in br[nb]) & set(partners):
                        shared = True
            rows.append({"domain": dom, "struct_pos": i,
                         "n_partners": len(partners),
                         "bridge_type": ("anti" if types == {"anti"} else
                                         "par" if types == {"par"} else
                                         "mixed" if types else "none"),
                         "has_bulge": int(shared)})
        if n % 200 == 0:
            print(f"  {n}/{len(paths)}", flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(out_csv, index=False)
    print(json.dumps({"strand_positions": int(len(df)),
                      "domains": int(df["domain"].nunique()),
                      "n_partners": {str(k): int(v) for k, v in
                                     df["n_partners"].value_counts().items()},
                      "bridge_type": {str(k): int(v) for k, v in
                                      df["bridge_type"].value_counts().items()},
                      "bulge_positions": int(df["has_bulge"].sum())}, indent=2))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
