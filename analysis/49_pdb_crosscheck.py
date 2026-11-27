"""A third, fully independent check on the secondary-structure assignment.

The assignment used throughout this project is a from-scratch implementation of the
Kabsch-Sander hydrogen-bond rules (analysis/17_dssp.py), run on AlphaFold models. It has
already been checked two ways: against backbone geometry the algorithm never uses
(18_ss_check.py), and against an independently produced assignment for the same domains
(23_ss_crosscheck.py, three-state agreement 0.801).

Both of those share something with the thing being tested: they are computed from the same
AlphaFold coordinates. This script removes that shared dependency. The Tsuboyama domains
are named by PDB identifier, so each has a deposited EXPERIMENTAL structure whose
HELIX and SHEET records were produced by the depositors, from different coordinates,
years before any of this. Comparing against those tests the implementation and the use of
predicted models at the same time.

Method. For each sampled domain: download the experimental entry from RCSB; read its CA
atoms per chain to get a sequence with author residue numbers; mark residues falling in a
HELIX or SHEET record. Separately read the AlphaFold model and run our assignment. The
domain sequence is then located inside a deposited chain by exact substring match -- no
alignment heuristics, so a match is unambiguous or the domain is skipped and counted.
Positions where the experimental structure has no CA (disordered, unresolved) are excluded,
because the depositor could not annotate what they could not see.

Reported: three-state agreement, the confusion matrix, and separately the helix/strand
confusion rate, which is the error that would actually matter -- loop/helix boundary
disagreement of a residue or two is expected between any two assignments and between a
crystal and a prediction.

Usage: python 49_pdb_crosscheck.py <afdir> <cache_dir> <out_json> [n_domains]
"""
import collections
import importlib.util
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("d17", os.path.join(HERE, "17_dssp.py"))
D17 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(D17)

THREE = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
         "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
         "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
         "TYR": "Y", "VAL": "V", "MSE": "M", "SEC": "U"}


def fetch(pdb_id, cache):
    path = os.path.join(cache, f"{pdb_id}.pdb")
    if os.path.exists(path) and os.path.getsize(path) > 1000:
        return path
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    r = subprocess.run(["curl", "-sL", "-m", "60", "-o", path, url], capture_output=True)
    if r.returncode != 0 or not os.path.exists(path) or os.path.getsize(path) < 1000:
        return None
    return path


def parse_experimental(path):
    """Return {chain: (sequence, [resnums], {resnum: 'H'|'E'})} from a deposited entry."""
    chains = collections.defaultdict(lambda: {"seq": [], "num": []})
    ss = collections.defaultdict(dict)
    seen = set()
    with open(path, errors="replace") as fh:
        for ln in fh:
            rec = ln[:6]
            if rec == "ENDMDL":
                break                      # first model only, for NMR entries
            if rec in ("ATOM  ", "HETATM") and ln[12:16].strip() == "CA":
                alt = ln[16]
                if alt not in (" ", "A"):
                    continue
                aa = THREE.get(ln[17:20].strip().upper())
                if aa is None:
                    continue
                ch = ln[21]
                key = (ch, ln[22:27])      # resSeq + icode
                if key in seen:
                    continue
                seen.add(key)
                try:
                    num = int(ln[22:26])
                except ValueError:
                    continue
                chains[ch]["seq"].append(aa)
                chains[ch]["num"].append(num)
            elif rec == "HELIX ":
                try:
                    ch, a, b = ln[19], int(ln[21:25]), int(ln[33:37])
                except ValueError:
                    continue
                for i in range(a, b + 1):
                    ss[ch][i] = "H"
            elif rec == "SHEET ":
                try:
                    ch, a, b = ln[21], int(ln[22:26]), int(ln[33:37])
                except ValueError:
                    continue
                for i in range(a, b + 1):
                    ss[ch][i] = "E"
    return {ch: ("".join(v["seq"]), v["num"], ss.get(ch, {})) for ch, v in chains.items()}


def main(afdir, cache, out_path, n_dom=200):
    os.makedirs(cache, exist_ok=True)
    models = sorted(f for f in os.listdir(afdir) if f.endswith(".pdb"))
    step = max(1, len(models) // int(n_dom))
    sample = models[::step][:int(n_dom)]

    conf = collections.Counter()
    skipped = collections.Counter()
    per_domain = []
    rows = []

    for fn in sample:
        pdb_id = fn[:-4].upper()
        model = os.path.join(afdir, fn)
        try:
            res = D17.parse_backbone(model)
            s, keys = D17.assign(res)
            if not s:
                skipped["model_unparsed"] += 1
                continue
            el = D17.elements(s, keys)
        except Exception:  # noqa: BLE001
            skipped["model_unparsed"] += 1
            continue
        ours_seq = "".join(res[i]["aa"] for i in keys)
        ours_ss = [el[i][0] for i in keys]

        exp_path = fetch(pdb_id, cache)
        if exp_path is None:
            skipped["not_downloadable"] += 1
            continue
        try:
            chains = parse_experimental(exp_path)
        except Exception:  # noqa: BLE001
            skipped["entry_unparsed"] += 1
            continue
        if not any(v[2] for v in chains.values()):
            skipped["no_ss_records"] += 1
            continue

        hit = None
        for ch, (seq, nums, ssmap) in chains.items():
            j = seq.find(ours_seq)
            if j >= 0:
                hit = (ch, nums[j:j + len(ours_seq)], ssmap)
                break
        if hit is None:
            skipped["sequence_not_found_in_entry"] += 1
            continue

        ch, nums, ssmap = hit
        agree = tot = 0
        for key, ours, num in zip(keys, ours_ss, nums):
            theirs = ssmap.get(num, "L")
            conf[(ours, theirs)] += 1
            tot += 1
            agree += (ours == theirs)
            rows.append((fn, key, ours, theirs))
        per_domain.append({"pdb": pdb_id, "chain": ch, "n": tot,
                           "agreement": round(agree / tot, 4) if tot else None})

    total = sum(conf.values())
    diag = sum(v for (a, b), v in conf.items() if a == b)
    he = sum(v for (a, b), v in conf.items() if {a, b} == {"H", "E"})
    rep = {
        "domains_sampled": len(sample),
        "domains_compared": len(per_domain),
        "skipped": dict(skipped),
        "residues_compared": total,
        "three_state_agreement": round(diag / total, 4) if total else None,
        "helix_strand_confusions": he,
        "helix_strand_confusion_rate": round(he / total, 5) if total else None,
        "confusion_matrix_ours_x_deposited": {f"{a}->{b}": v for (a, b), v in
                                              sorted(conf.items())},
        "per_class_recall": {
            c: round(conf[(c, c)] / max(sum(v for (a, b), v in conf.items() if a == c), 1), 4)
            for c in "HEL"},
        "note": "deposited HELIX/SHEET records come from the experimental structure and "
                "are independent of both this implementation and the AlphaFold models; "
                "residues absent from the records are treated as loop, which is how the "
                "format expresses them",
        "median_per_domain_agreement": (
            round(sorted(d["agreement"] for d in per_domain)[len(per_domain) // 2], 4)
            if per_domain else None)}
    json.dump(rep, open(out_path, "w"), indent=2)
    # per-residue table, so downstream scripts can restrict to positions where
    # our assignment and the deposited structure agree
    with open(out_path.replace(".json", "_residues.csv"), "w", newline="") as fh:
        fh.write("domain,pos,ours,deposited\n")
        for d, k, a, b in rows:
            fh.write(f"{d},{k},{a},{b}\n")

    print(f"domains compared: {rep['domains_compared']} / {rep['domains_sampled']}")
    print(f"skipped: {rep['skipped']}")
    print(f"residues: {rep['residues_compared']}")
    print(f"three-state agreement: {rep['three_state_agreement']}")
    print(f"median per-domain agreement: {rep['median_per_domain_agreement']}")
    print(f"helix/strand confusions: {rep['helix_strand_confusions']} "
          f"({rep['helix_strand_confusion_rate']})")
    print("confusion (ours -> deposited):")
    for k, v in rep["confusion_matrix_ours_x_deposited"].items():
        print(f"   {k}  {v}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3],
         sys.argv[4] if len(sys.argv) > 4 else 200)
