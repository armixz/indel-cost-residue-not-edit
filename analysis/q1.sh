#!/usr/bin/env bash
PY=$(dirname "$(readlink -f /root/.local/bin/state)")/python
exec "$PY" -I - "$@" <<'PY'
import sys, re, numpy as np, pandas as pd
from scipy.stats import spearmanr
df = pd.read_csv("/mnt/c/Users/armin/Desktop/AXZ/insertion-lab/data/mavedb/kir21_scores.csv")
rows=[]
for p,s in zip(df["hgvs_pro"].fillna(""), df["score"]):
    if not np.isfinite(s): continue
    m=re.match(r"p\.([A-Za-z]{3})(\d+)_[A-Za-z]{3}\d+ins([A-Za-z]+)$",p)
    if m:
        n=len(re.findall(r"[A-Z][a-z]{2}",m.group(3)))
        rows.append((int(m.group(2)),f"ins{n}",float(s)))
    d=re.match(r"p\.[A-Za-z]{3}(\d+)del$",p)
    if d: rows.append((int(d.group(1)),"del1",float(s)))
t=pd.DataFrame(rows,columns=["pos","kind","score"]).pivot_table(index="pos",columns="kind",values="score")
rng=np.random.default_rng(0)
def ci(x,y,n=2000):
    v=[spearmanr(x[i],y[i]).statistic for i in (rng.integers(0,len(x),len(x)) for _ in range(n))]
    return [round(float(q),3) for q in np.nanquantile(v,[.025,.975])]
print("positions:",len(t),"coverage:",{c:int(t[c].notna().sum()) for c in t.columns})
print("medians:",{c:round(float(t[c].median()),3) for c in t.columns})
for a,b in [("ins1","ins2"),("ins2","ins3"),("ins1","ins3"),("ins1","del1"),("ins3","del1")]:
    m=t[a].notna()&t[b].notna(); x,y=t.loc[m,a].values,t.loc[m,b].values
    print(f"{a} vs {b}: n={m.sum()} spearman={spearmanr(x,y).statistic:.3f} ci95={ci(x,y)}")
# does the SET of tolerated sites change with size?
for q in (0.25,0.5):
    thr1,thr3=t["ins1"].quantile(q),t["ins3"].quantile(q)
    m=t["ins1"].notna()&t["ins3"].notna()
    bad1,bad3=t.loc[m,"ins1"]<thr1,t.loc[m,"ins3"]<thr3
    inter=(bad1&bad3).sum(); union=(bad1|bad3).sum()
    print(f"bottom {int(q*100)}% by own threshold: jaccard overlap 1aa vs 3aa = {inter/union:.3f} ({inter}/{union})")
PY
