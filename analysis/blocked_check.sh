#!/usr/bin/env bash
PY=$(dirname "$(readlink -f /root/.local/bin/state)")/python
exec "$PY" -I - <<'PY'
import importlib.util, numpy as np, json
from scipy.stats import rankdata
from sklearn.metrics import roc_auc_score
root="/mnt/c/Users/armin/Desktop/AXZ/insertion-lab"
spec=importlib.util.spec_from_file_location("g6",f"{root}/analysis/10_gate6_how_many.py")
g6=importlib.util.module_from_spec(spec); spec.loader.exec_module(g6)
data=g6.build_datasets(root)
ORI={"plddt":-1,"rsa":+1,"dist_term":-1,"prodomino":+1}
CAND=["plddt","rsa","dist_term","prodomino","consensus"]
def auroc(y,s):
    m=np.isfinite(s)
    if m.sum()<5 or len(np.unique(y[m]))<2: return np.nan
    return float(roc_auc_score(y[m],s[m]))
def prep(d):
    cols={}
    for f in CAND:
        if f=="consensus":
            base=[c for c in ("plddt","rsa","dist_term","prodomino") if c in d]
            cols[f]=np.mean([rankdata(ORI[c]*np.nan_to_num(d[c].values.astype(float),nan=float(np.nanmedian(d[c])))) for c in base],0)
        else:
            cols[f]=ORI[f]*np.nan_to_num(d[f].values.astype(float),nan=float(np.nanmedian(d[f])))
    return cols
K=24; REP=50
res={}
for nm,d in data.items():
    y=d["y"].values.astype(int); cand=prep(d); n=len(y)
    order=np.argsort(d["pos"].values)          # contiguous along the sequence
    sel=[];fx={c:[] for c in CAND};orc=[]
    for r in range(REP):
        rng=np.random.default_rng(99000+r)
        # panel = 3 contiguous blocks of 8 consecutive positions
        tr=[]
        for _ in range(3):
            st=int(rng.integers(0,max(1,n-8)))
            tr+= list(order[st:st+8])
        tr=np.unique(tr); te=np.setdiff1d(np.arange(n),tr)
        if len(np.unique(y[tr]))<2 or len(np.unique(y[te]))<2: continue
        panel={c:auroc(y[tr],cand[c][tr]) for c in CAND}
        fin=[(v,c) for c,v in panel.items() if np.isfinite(v)]
        if not fin: continue
        best=max(fin)[1]
        sel.append(auroc(y[te],cand[best][te]))
        for c in CAND: fx[c].append(auroc(y[te],cand[c][te]))
        orc.append(max(auroc(y[te],cand[c][te]) for c in CAND))
    if sel:
        res[nm]={"select":round(float(np.nanmedian(sel)),3),
                 "plddt":round(float(np.nanmedian(fx["plddt"])),3),
                 "prodomino":round(float(np.nanmedian(fx["prodomino"])),3),
                 "oracle":round(float(np.nanmedian(orc)),3),"n_rep":len(sel)}
print("BLOCKED panel (3 contiguous runs of 8), k=24")
for nm,v in res.items(): print(f"  {nm:16s} {v}")
for m in ("select","plddt","prodomino","oracle"):
    v=[res[n][m] for n in res]
    print(f"{m:10s} mean={np.mean(v):.3f} worst={min(v):.3f}")
PY
