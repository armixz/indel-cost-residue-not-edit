"""Gate 0: how much measured domain-insertion data is actually usable?

Counts, per protein and inserted domain, how many distinct in-frame insertion sites have reads in
both the unsorted input library and the sorted gates, which is what a per-site enrichment score
needs. Reads only; computes no model and fits nothing.

Usage: python -I 00_census.py <path to DI_screen/analysis/input_data> <out.json>
"""
import json
import os
import re
import sys
from collections import defaultdict

import pandas as pd

COLS = ["insertion_site", "in_frame_insertion", "forward_insertion", "valid_linker"]


def load(path):
    df = pd.read_csv(path, usecols=lambda c: c in COLS, low_memory=False)
    for c in ("in_frame_insertion", "forward_insertion", "valid_linker"):
        if c in df:
            df = df[df[c].astype(str).str.lower().isin(("true", "1", "1.0"))]
    return df["insertion_site"].dropna().astype(int)


def main(src, out_path):
    files = sorted(f for f in os.listdir(src) if f.endswith("_combined.csv"))
    groups = defaultdict(dict)
    for f in files:
        m = re.match(r"(?P<prot>[^_]+)_(?P<dom>[^_]+)_(?P<gate>input|12|22|22_ind)_combined\.csv", f)
        if not m:
            print("unparsed:", f, file=sys.stderr)
            continue
        groups[(m["prot"], m["dom"])][m["gate"]] = os.path.join(src, f)

    report = {"source": src, "groups": {}}
    for (prot, dom), gates in sorted(groups.items()):
        entry = {"gates": {}}
        sites = {}
        for gate, path in sorted(gates.items()):
            s = load(path)
            vc = s.value_counts()
            sites[gate] = vc
            entry["gates"][gate] = {
                "reads_in_frame": int(len(s)),
                "distinct_sites": int(vc.size),
                "sites_ge10_reads": int((vc >= 10).sum()),
                "site_min": int(s.min()) if len(s) else None,
                "site_max": int(s.max()) if len(s) else None,
                "file_mb": round(os.path.getsize(path) / 1048576, 1),
            }
        if "input" in sites:
            for gate in [g for g in sites if g != "input"]:
                both = sites["input"].index.intersection(sites[gate].index)
                deep = sites["input"][sites["input"] >= 10].index.intersection(
                    sites[gate][sites[gate] >= 10].index)
                entry.setdefault("usable_sites", {})[gate] = {
                    "in_input_and_gate": int(len(both)),
                    "both_ge10_reads": int(len(deep)),
                }
        report["groups"][f"{prot}_{dom}"] = entry

    prots = sorted({k.split("_")[0] for k in report["groups"]})
    best = {}
    for k, v in report["groups"].items():
        p = k.split("_")[0]
        for gate, u in (v.get("usable_sites") or {}).items():
            best[p] = max(best.get(p, 0), u["both_ge10_reads"])
    report["summary"] = {
        "proteins": prots, "n_proteins": len(prots),
        "protein_domain_pairs": len(report["groups"]),
        "best_usable_sites_per_protein_ge10_reads": best,
        "total_usable_sites_ge10_reads": int(sum(best.values())),
    }
    json.dump(report, open(out_path, "w"), indent=1)
    print(json.dumps(report["summary"], indent=1))
    for k, v in report["groups"].items():
        print(k, {g: (d["distinct_sites"], d["sites_ge10_reads"]) for g, d in v["gates"].items()},
              "usable:", v.get("usable_sites"))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
