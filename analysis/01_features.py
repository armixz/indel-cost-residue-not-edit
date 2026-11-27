"""Labels and baseline features for the insertion-tolerance benchmark (see ../RESULTS.md).

Label follows DI_screen/analysis/1-sequencing_processing.ipynb exactly.
Features: conservation entropy (provided MSA), relative SASA and secondary structure and pLDDT
(provided AlphaFold model), distance to the nearest terminus.

Usage: python -I 01_features.py <DI_screen/analysis/input_data> <out.csv>
"""
import math
import os
import sys

import numpy as np
import pandas as pd

SAMPLES = ["AraC_ERD", "AraC_eYFP", "AraC_LOV", "AraC_PDZ", "AraC_UniRapR", "Flp_PDZ", "SigF_PDZ", "TVMV_PDZ"]
REPS = ("12", "22")
# Tien et al. 2013 theoretical maximum accessible surface area, A^2
MAXASA = {"A": 129, "R": 274, "N": 195, "D": 193, "C": 167, "E": 223, "Q": 225, "G": 104, "H": 224,
          "I": 197, "L": 201, "K": 236, "M": 224, "F": 240, "P": 159, "S": 155, "T": 172, "W": 285,
          "Y": 263, "V": 174}
VDW = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80}
THREE2ONE = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLU": "E", "GLN": "Q",
             "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
             "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}


def read_fasta(path):
    seqs, name = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = []
        elif name:
            seqs[name].append(line)
    return {k: "".join(v) for k, v in seqs.items()}


# ---------- label, following the authors' notebook ----------

def site_counts(path):
    df = pd.read_csv(path, low_memory=True,
                     usecols=["insertion_site", "in_frame_insertion", "forward_insertion", "linker_seq"])
    keep = (df["in_frame_insertion"].astype(str).str.lower() == "true") & \
           (df["forward_insertion"].astype(str).str.lower() == "true") & df["linker_seq"].isna()
    s = df.loc[keep, "insertion_site"] / 3.0 - 1.0
    s = s[np.isfinite(s)]
    return s.astype(int).value_counts()


def enrichment(src, sample, length):
    """Mean over the two replicates of log2(share in gate / share in input), authors' rules."""
    inp = site_counts(os.path.join(src, f"{sample}_input_combined.csv"))
    inp = inp / inp.sum()
    out = {}
    for rep in REPS:
        p = os.path.join(src, f"{sample}_{rep}_combined.csv")
        if not os.path.exists(p):
            continue
        c = site_counts(p)
        c = c / c.sum()
        idx = pd.RangeIndex(0, length)
        cond = c.reindex(idx).fillna(0.0)
        base = inp.reindex(idx)
        log = np.log2(cond / base)
        log[~np.isfinite(log) & (cond == 0)] = -10.0      # authors set -inf to -10
        log[base.isna() | (base == 0)] = np.nan           # absent from input -> missing
        out[rep] = log
    if not out:
        return None
    reps = pd.DataFrame(out)
    depth_in = inp.reindex(pd.RangeIndex(0, length)).fillna(0.0)
    raw_in = site_counts(os.path.join(src, f"{sample}_input_combined.csv")).reindex(pd.RangeIndex(0, length)).fillna(0)
    return pd.DataFrame({"score": reps.mean(axis=1, skipna=False), "n_reps": reps.notna().sum(axis=1),
                         "input_share": depth_in, "input_reads": raw_in})


# ---------- structure features ----------

def parse_pdb(path):
    res = {}
    for line in open(path):
        if not line.startswith("ATOM"):
            continue
        name = line[12:16].strip()
        rn = THREE2ONE.get(line[17:20].strip())
        if rn is None or line[16] not in (" ", "A"):
            continue
        i = int(line[22:26])
        el = (line[76:78].strip() or name[0]).upper()
        xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        b = float(line[60:66])
        res.setdefault(i, {"aa": rn, "atoms": [], "b": []})
        res[i]["atoms"].append((name, el, xyz))
        res[i]["b"].append(b)
    return res


def sphere_points(n=92):
    pts, off, inc = [], 2.0 / n, math.pi * (3.0 - math.sqrt(5.0))
    for k in range(n):
        y = k * off - 1 + off / 2
        r = math.sqrt(max(0.0, 1 - y * y))
        phi = k * inc
        pts.append((math.cos(phi) * r, y, math.sin(phi) * r))
    return np.array(pts)


def shrake_rupley(res, probe=1.4, n=92):
    atoms, owner = [], []
    for i, d in sorted(res.items()):
        for nm, el, xyz in d["atoms"]:
            atoms.append((xyz, VDW.get(el, 1.7) + probe))
            owner.append(i)
    xyz = np.array([a[0] for a in atoms])
    rad = np.array([a[1] for a in atoms])
    owner = np.array(owner)
    pts = sphere_points(n)
    sasa = np.zeros(len(atoms))
    for k in range(len(atoms)):
        d = np.linalg.norm(xyz - xyz[k], axis=1)
        nb = np.flatnonzero((d < rad + rad[k]) & (d > 0))
        test = xyz[k] + rad[k] * pts
        if len(nb):
            dist = np.linalg.norm(test[:, None, :] - xyz[nb][None, :, :], axis=2)
            acc = np.all(dist >= rad[nb][None, :], axis=1)
        else:
            acc = np.ones(len(test), bool)
        sasa[k] = 4 * math.pi * rad[k] ** 2 * acc.mean()
    out = {}
    for i in np.unique(owner):
        out[int(i)] = float(sasa[owner == i].sum())
    return out


def dihedral(p0, p1, p2, p3):
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1n = b1 / (np.linalg.norm(b1) + 1e-9)
    v = b0 - np.dot(b0, b1n) * b1n
    w = b2 - np.dot(b2, b1n) * b1n
    return math.degrees(math.atan2(np.dot(np.cross(b1n, v), w), np.dot(v, w)))


def backbone_ss(res):
    """Helix / strand / loop from phi and psi (DSSP is not available here)."""
    get = lambda i, nm: next((np.array(x[2]) for x in res[i]["atoms"] if x[0] == nm), None) if i in res else None
    ss = {}
    for i in sorted(res):
        C0, N1, CA1, C1, N2 = get(i - 1, "C"), get(i, "N"), get(i, "CA"), get(i, "C"), get(i + 1, "N")
        if any(v is None for v in (C0, N1, CA1, C1, N2)):
            ss[i] = "L"
            continue
        phi, psi = dihedral(C0, N1, CA1, C1), dihedral(N1, CA1, C1, N2)
        if -160 <= phi <= -20 and -120 <= psi <= 50:
            ss[i] = "H"
        elif -180 <= phi <= -40 and (90 <= psi <= 180 or -180 <= psi <= -150):
            ss[i] = "E"
        else:
            ss[i] = "L"
    return ss


def msa_conservation(path, ref_seq):
    seqs = read_fasta(path)
    names = list(seqs)
    ref = next((n for n in names if seqs[n].replace("-", "").upper() == ref_seq.upper()), None)
    if ref is None:                                   # fall back to the best-matching sequence
        def score(n):
            u = seqs[n].replace("-", "").upper()
            return sum(a == b for a, b in zip(u, ref_seq.upper())) / max(len(ref_seq), 1)
        ref = max(names, key=score)
    aln = np.array([list(seqs[n].upper()) for n in names])
    col_of_res, r = {}, 0
    for c, ch in enumerate(seqs[ref].upper()):
        if ch != "-":
            col_of_res[r] = c
            r += 1
    out = {}
    for r, c in col_of_res.items():
        col = aln[:, c]
        aas = col[(col != "-") & np.isin(col, list(MAXASA))]
        if len(aas) < 2:
            out[r] = (np.nan, float((col == "-").mean()))
            continue
        _, cnt = np.unique(aas, return_counts=True)
        p = cnt / cnt.sum()
        out[r] = (float(-(p * np.log2(p)).sum()), float((col == "-").mean()))
    return out, ref, len(names)


def main(src, out_path):
    prots = read_fasta(os.path.join(src, "proteins.fasta"))
    rows = []
    cache = {}
    for sample in SAMPLES:
        prot = sample.split("_")[0]
        seq = prots[prot]
        L = len(seq)
        if prot not in cache:
            res = parse_pdb(os.path.join(src, f"{prot}_AF.pdb"))
            sasa = shrake_rupley(res)
            ss = backbone_ss(res)
            cons, refname, nseq = msa_conservation(os.path.join(src, "MSA", f"{prot}_ref_aln.afa"), seq)
            plddt = {i: float(np.mean(d["b"])) for i, d in res.items()}
            cache[prot] = (res, sasa, ss, cons, plddt, refname, nseq)
            print(f"{prot}: L={L} pdb_res={len(res)} msa_seqs={nseq} ref={refname}", flush=True)
        res, sasa, ss, cons, plddt, refname, nseq = cache[prot]
        e = enrichment(src, sample, L)
        if e is None:
            continue
        for r in range(L):
            i = r + 1                                   # PDB numbering is 1-based
            ent, gap = cons.get(r, (np.nan, np.nan))
            aa = res.get(i, {}).get("aa")
            rsa = sasa.get(i)
            rows.append({
                "protein": prot, "domain": sample.split("_", 1)[1], "res0": r, "aa": aa,
                "score": e["score"].get(r, np.nan), "n_reps": int(e["n_reps"].get(r, 0)),
                "input_reads": int(e["input_reads"].get(r, 0)),
                "entropy": ent, "gap_frac": gap,
                "rsa": (rsa / MAXASA[aa]) if (rsa is not None and aa in MAXASA) else np.nan,
                "ss": ss.get(i), "plddt": plddt.get(i),
                "dist_term": min(r, L - 1 - r), "length": L,
            })
    df = pd.DataFrame(rows)
    df.to_csv(out_path, index=False)
    print("\nwrote", out_path, df.shape)
    g = df.dropna(subset=["score"]).groupby(["protein", "domain"])
    print(g.agg(sites=("score", "size"), positive=("score", lambda s: int((s > 0).sum())),
                mean_score=("score", "mean")).to_string())


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
