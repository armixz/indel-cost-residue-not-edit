#!/usr/bin/env bash
PY=$(dirname "$(readlink -f /root/.local/bin/state)")/python
exec "$PY" -I - <<'PY'
import pandas as pd, numpy as np
d=pd.read_csv("/mnt/c/Users/armin/Desktop/AXZ/insertion-lab/results/gate11_positions.csv")
print("helix vs loop delta WITHIN each burial tertile:")
for b,g in d.groupby("burial"):
    h=g[g.ss=="H"]["delta"]; l=g[g.ss=="L"]["delta"]
    if len(h)>30 and len(l)>30:
        print(f"  {b:8s} helix {h.mean():+.3f} (n={len(h)})  loop {l.mean():+.3f} (n={len(l)})  diff {h.mean()-l.mean():+.3f}")
print("\nburial effect WITHIN each ss class:")
for s,g in d.groupby("ss"):
    bu=g[g.burial=="buried"]["delta"]; ex=g[g.burial=="exposed"]["delta"]
    if len(bu)>30 and len(ex)>30:
        print(f"  {s} buried {bu.mean():+.3f} exposed {ex.mean():+.3f} diff {bu.mean()-ex.mean():+.3f}")
print("\npositive-phi flank effect WITHIN each ss class:")
for s,g in d.groupby("ss"):
    a=g[g.flank_pos_phi==1]["delta"]; b=g[g.flank_pos_phi==0]["delta"]
    if len(a)>30 and len(b)>30:
        print(f"  {s} flankPosPhi {a.mean():+.3f} (n={len(a)})  other {b.mean():+.3f}  diff {a.mean()-b.mean():+.3f}")
# post-hoc contrast: helix middle vs loop
hm=d[(d.ss=="H")&(d.helix_middle==1)]["delta"]; lo=d[d.ss=="L"]["delta"]
print(f"\npost-hoc: helix-middle {hm.mean():+.3f} vs loop {lo.mean():+.3f} -> {hm.mean()-lo.mean():+.3f} kcal/mol")
# domain-level bootstrap for that contrast
rng=np.random.default_rng(0); doms=d.domain.unique()
idx={x:d.index[d.domain==x].to_numpy() for x in doms}
v=[]
for _ in range(2000):
    s=d.loc[np.concatenate([idx[x] for x in rng.choice(doms,len(doms),True)])]
    a=s[(s.ss=="H")&(s.helix_middle==1)]["delta"].mean(); b=s[s.ss=="L"]["delta"].mean()
    v.append(a-b)
print("  ci95 by domain:", [round(float(q),3) for q in np.nanquantile(v,[.025,.975])])
PY
