"""Extract single-residue insertions and deletions from the Tsuboyama et al. 2023 Lib3 table.

Names look like `1E0L.pdb_insG12` (insert Gly at position 12), `1E0L.pdb_insA12` (Ala),
`1E0L.pdb_delT3` (delete the Thr at position 3). Only wild-type backgrounds are kept: anything
with an extra substitution before the indel token is dropped.

Usage: python -I 13_tsuboyama_extract.py <Lib3_K50dG.csv> <out.csv>
"""
import csv
import re
import sys

PAT_INS = re.compile(r"^([0-9A-Za-z_]+)\.pdb_ins([GA])(\d+)$")
PAT_DEL = re.compile(r"^([0-9A-Za-z_]+)\.pdb_del([A-Z])(\d+)$")
PAT_WT = re.compile(r"^([0-9A-Za-z]{4})$")


def main(src, out):
    n = kept = 0
    with open(src, newline="") as fh, open(out, "w", newline="") as oh:
        rd = csv.DictReader(fh)
        wr = csv.writer(oh)
        wr.writerow(["domain", "kind", "res", "pos", "dG", "dG_lo", "dG_hi", "dG_ci"])
        for row in rd:
            n += 1
            name = row["name"]
            m = PAT_INS.match(name)
            kind = "ins"
            if not m:
                m = PAT_DEL.match(name)
                kind = "del"
            if not m:
                continue
            dom = m.group(1)
            if "_" in dom:                       # a mutant background, not wild type
                continue
            try:
                dg = float(row["deltaG"]); lo = float(row["deltaG_95CI_low"])
                hi = float(row["deltaG_95CI_high"]); ci = float(row["deltaG_95CI"])
            except (ValueError, KeyError, TypeError):
                continue
            wr.writerow([dom, kind, m.group(2), int(m.group(3)), dg, lo, hi, ci])
            kept += 1
            if kept % 20000 == 0:
                print("kept", kept, "of", n, flush=True)
    print("rows read", n, "kept", kept)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
