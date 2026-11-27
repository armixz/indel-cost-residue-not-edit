"""Could H-strand be inherited from the unfolded-state model, or from the assay itself?

Tsuboyama et al. infer dG from protease resistance: dG is a function of the measured K50 and of a
PREDICTED unfolded-state K50 computed from sequence. A beta strand's sequence neighbourhood is
hydrophobic and a helix's is not, so if the unfolded-state model is systematically off by
structure class, the dG difference reported here would be an artefact of that model rather than a
property of the folded protein.

Checks:
  1. the paired strand-minus-helix difference in log10_K50unfolded itself (trypsin and
     chymotrypsin). This is pure model output with no folded-state information in it, so under the
     artefact explanation it should show the same contrast.
  2. the same paired difference in the raw measured log10_K50, which contains no unfolded-state
     model at all. If H-strand is real it should appear here too.
  3. consistency between the two proteases: an assay artefact would not have to agree across them.
  4. the three domains shared with the abundancePCA data, assay against assay at matched positions.

Usage: indel_env python 32_unfolded_check.py <lib3_csv> <results_dir> <apca_csv> <out_json>
"""
import csv
import json
import sys

import numpy as np
import pandas as pd
from scipy import stats

RNG = np.random.default_rng(0)
SHARED = {"CSPA-CSD": "1MJC", "CSPB-CSD": "1CSQ", "FBP11-FF1": "1UZC"}


def agg(vals, nboot=20000):
    v = np.asarray([x for x in vals if x is not None and np.isfinite(x)], dtype=float)
    if len(v) < 3:
        return {"n_domains": int(len(v)), "mean": None, "ci95": [None, None]}
    idx = RNG.integers(0, len(v), size=(nboot, len(v)))
    bs = v[idx].mean(axis=1)
    return {"n_domains": int(len(v)), "mean": round(float(v.mean()), 3),
            "ci95": [round(float(np.quantile(bs, 0.025)), 3),
                     round(float(np.quantile(bs, 0.975)), 3)]}


def paired_EH(df, col, nmin=2):
    out = []
    for _, g in df.groupby("domain"):
        e, h = g[g["ss"] == "E"][col], g[g["ss"] == "H"][col]
        if len(e) >= nmin and len(h) >= nmin:
            out.append(float(e.mean() - h.mean()))
    return out


def read_lib3(path):
    rows = []
    csv.field_size_limit(10 ** 7)
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh):
            nm = r["name"]
            if "_ins" not in nm and "_del" not in nm:
                continue
            kind = "ins" if "_ins" in nm else "del"
            sep = ".pdb_" if ".pdb_" in nm else "_"
            parts = nm.split(sep, 1)
            if len(parts) != 2:
                continue
            dom, tail = parts
            i = tail.find("ins" if kind == "ins" else "del")
            if i < 0:
                continue
            tail = tail[i + 3:]
            res, num = tail[0], tail[1:]
            if not num.isdigit():
                continue
            def f(k):
                try:
                    return float(r[k])
                except (ValueError, KeyError):
                    return np.nan
            rows.append({"domain": dom, "kind": kind, "res": res, "pos": int(num),
                         "dG": f("deltaG"), "K50_t": f("log10_K50_t"), "K50_c": f("log10_K50_c"),
                         "U_t": f("log10_K50unfolded_t"), "U_c": f("log10_K50unfolded_c")})
    return pd.DataFrame(rows)


def main(lib3, rdir, apca_csv, out_path):
    raw = read_lib3(lib3)
    feat = pd.read_csv(f"{rdir}/tsub_features.csv")
    ss = pd.read_csv(f"{rdir}/tsub_ss.csv")
    ss["domain"] = ss["domain"].str.replace(".pdb", "", regex=False)
    base = feat[["domain", "pos", "struct_pos"]].merge(
        ss[["domain", "struct_pos", "ss"]], on=["domain", "struct_pos"])

    ins = raw[(raw["kind"] == "ins") & (raw["res"] == "A")]
    dele = raw[raw["kind"] == "del"]
    rep = {"n_raw_rows": int(len(raw))}

    for label, sub in (("insertion_Ala", ins), ("deletion", dele)):
        m = sub.merge(base, on=["domain", "pos"], how="inner")
        m = m[m["ss"].isin(["E", "H"])]
        out = {"n_positions": int(len(m)),
               "n_domains_with_both": len(paired_EH(m, "dG"))}
        for col, nice in (("dG", "deltaG"), ("K50_t", "raw_K50_trypsin"),
                          ("K50_c", "raw_K50_chymotrypsin"),
                          ("U_t", "UNFOLDED_model_trypsin"),
                          ("U_c", "UNFOLDED_model_chymotrypsin")):
            mm = m[m[col].notna()]
            out[nice] = agg(paired_EH(mm, col))
        rep[label] = out

    # 3. protease consistency on dG is implicit; check correlation of the two K50 contrasts
    m = ins.merge(base, on=["domain", "pos"], how="inner")
    m = m[m["ss"].isin(["E", "H"]) & m["K50_t"].notna() & m["K50_c"].notna()]
    vt, vc = paired_EH(m, "K50_t"), paired_EH(m, "K50_c")
    n = min(len(vt), len(vc))
    rep["protease_agreement"] = {
        "n_domains": n,
        "pearson_trypsin_vs_chymotrypsin_contrast":
            round(float(stats.pearsonr(vt[:n], vc[:n])[0]), 3) if n > 5 else None}

    # 4. shared domains, assay vs assay
    ap = pd.read_csv(apca_csv)
    ap = ap[ap["mut_type"].isin(["insertions", "deletions"])].copy()
    ap["cost_apca"] = -ap["scaled_fitness"]
    ap["ss_apca"] = np.where(ap["Structure_name"].isin(["Strand", "Bridge"]), "E",
                             np.where(ap["Structure_name"].isin(["AlphaHelix", "310Helix"]), "H", "L"))
    shared = {}
    for aname, tname in SHARED.items():
        a = ap[ap["domain"] == aname]
        t = raw[raw["domain"] == tname]
        if len(a) == 0 or len(t) == 0:
            shared[aname] = {"status": "missing"}
            continue
        tb = base[base["domain"] == tname]
        ti = (t[t["kind"] == "ins"].merge(tb, on=["domain", "pos"])
              .groupby("struct_pos")["dG"].mean())
        ai = (a[a["mut_type"] == "insertions"].groupby("Pos")["cost_apca"].mean())
        # align by offset search: the two numbering systems may differ by a constant
        best = None
        for off in range(-10, 11):
            idx = [p for p in ti.index if (p + off) in ai.index]
            if len(idx) < 15:
                continue
            x = np.array([ti[p] for p in idx])
            y = np.array([ai[p + off] for p in idx])
            r = float(stats.spearmanr(x, -y)[0])  # dG high = stable; cost high = bad
            if best is None or r > best[1]:
                best = (off, r, len(idx))
        shared[aname] = {"tsuboyama_domain": tname,
                         "best_offset": best[0] if best else None,
                         "spearman_dG_vs_minus_apca_cost": round(best[1], 3) if best else None,
                         "n_matched_positions": best[2] if best else 0}
    rep["shared_domain_assay_agreement"] = shared

    json.dump(rep, open(out_path, "w"), indent=2)
    print(json.dumps(rep, indent=2))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
