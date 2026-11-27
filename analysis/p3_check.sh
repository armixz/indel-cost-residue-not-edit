#!/usr/bin/env bash
PY=$(dirname "$(readlink -f /root/.local/bin/state)")/python
exec "$PY" -I - <<'PY'
import sys, re, numpy as np, pandas as pd
from scipy.stats import spearmanr, rankdata
sys.path.insert(0,"/mnt/c/Users/armin/Desktop/AXZ/insertion-lab/analysis")
md="/mnt/c/Users/armin/Desktop/AXZ/insertion-lab/data/mavedb"
import importlib.util
spec=importlib.util.spec_from_file_location("g4","/mnt/c/Users/armin/Desktop/AXZ/insertion-lab/analysis/08_gate4_segments.py")
g4=importlib.util.module_from_spec(spec); spec.loader.exec_module(g4)
prof=g4.plddt_profile(f"{md}/af_P35561.pdb"); rl,rs,thr=g4.runs(prof,0.5)
a=g4.mavedb_insertions(f"{md}/kir21_scores.csv",1).set_index("pos")
b=g4.mavedb_insertions(f"{md}/kir21_scores.csv",3).set_index("pos")
common=[p for p in a.index.intersection(b.index) if g4.kir_map(int(p)) in prof]
st=[g4.kir_map(int(p)) for p in common]
sens=(a.loc[common,"score"]-b.loc[common,"score"]).values
slack=np.array([rs[s] for s in st],float); pl=np.array([prof[s] for s in st],float)
runl=np.array([rl[s] for s in st],float)
def pres(x,z):
    X=np.c_[np.ones(len(z)),rankdata(z)]; r=rankdata(x)
    return r-X@np.linalg.lstsq(X,r,rcond=None)[0]
print("n",len(common))
print("raw   slack vs size-sensitivity : %.3f"%spearmanr(slack,sens).statistic)
print("raw   pLDDT vs size-sensitivity : %.3f"%spearmanr(pl,sens).statistic)
print("raw   run_len vs size-sensitivity: %.3f"%spearmanr(runl,sens).statistic)
print("raw   slack vs pLDDT            : %.3f"%spearmanr(slack,pl).statistic)
pr=spearmanr(pres(slack,pl),pres(sens,pl)).statistic
rng=np.random.default_rng(0)
bs=[spearmanr(pres(slack[i],pl[i]),pres(sens[i],pl[i])).statistic for i in (rng.integers(0,len(sens),len(sens)) for _ in range(2000))]
print("PARTIAL slack vs sensitivity | pLDDT : %.3f  ci95=%s"%(pr,[round(float(q),3) for q in np.nanquantile(bs,[.025,.975])]))
# and the reverse: pLDDT controlling for slack
pr2=spearmanr(pres(pl,slack),pres(sens,slack)).statistic
print("PARTIAL pLDDT vs sensitivity | slack : %.3f"%pr2)
PY
