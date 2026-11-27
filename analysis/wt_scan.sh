#!/usr/bin/env bash
PY=$(dirname "$(readlink -f /root/.local/bin/state)")/python
cd "/mnt/c/Users/armin/Desktop/AXZ/insertion-lab/data/zenodo/tsub_dg/K50_dG_tables"
exec "$PY" -I - <<'PY'
import csv,re,json,statistics as st
R="/mnt/c/Users/armin/Desktop/AXZ/insertion-lab/results"
want={r["domain"] for r in csv.DictReader(open(f"{R}/tsuboyama_indels.csv"))}
wt={}
for lib in ("Lib1","Lib2","Lib3","Lib4"):
    with open(f"{lib}_K50dG.csv",newline="") as fh:
        for row in csv.DictReader(fh):
            m=re.match(r"^([0-9A-Za-z]{4})(?:\.pdb)?$",row["name"])
            if m and m.group(1) in want:
                try: wt.setdefault(m.group(1),[]).append(float(row["deltaG"]))
                except: pass
print("domains wanted",len(want),"with WT dG",len(wt))
print("with >1 measurement:",sum(1 for v in wt.values() if len(v)>1))
json.dump({k:st.median(v) for k,v in wt.items()},open(f"{R}/tsub_wt_dG.json","w"))
PY
