"""R1 rerun on the COMPLETE set of 748 domains.

The first pass used a loader that merged a cached secondary-structure table and silently
dropped 57 HEEH_KT designs whose names carry no .pdb suffix. This repeats the natural
versus designed comparison using the same loader the manuscript uses, so the designed
subset is complete (349 rather than 292).
"""
import importlib.util, json, os, re, sys
import numpy as np, pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
def mod(n, f):
    sp = importlib.util.spec_from_file_location(n, os.path.join(HERE, f))
    m = importlib.util.module_from_spec(sp); sp.loader.exec_module(m); return m
D43 = mod("d43", "43_discriminate.py"); R52 = mod("r52", "52_robustness.py")
PACE = dict(D43.PACE); RNG = np.random.default_rng(0)
NAT = re.compile(r"^[0-9][A-Za-z0-9]{3}(\.pdb)?$")

lib3, afdir, rdir, out = sys.argv[1:5]
d = D43.load(lib3, afdir, rdir)
d["kind"] = ["natural" if NAT.match(x) else "designed" for x in d["domain"]]

def cells(sub):
    c = {}
    for (dm, aa), g in sub.groupby(["domain", "aa"]):
        h, l = g[g["ss"] == "H"]["cost"], g[g["ss"] == "L"]["cost"]
        if len(h) and len(l):
            c.setdefault(str(dm), {})[str(aa)] = float(h.mean() - l.mean())
    return c

def rho_from(cell, doms, min_dom=5):
    acc = {}
    for dm in doms:
        for aa, v in cell.get(dm, {}).items():
            acc.setdefault(aa, []).append(v)
    ctr = {a: float(np.mean(v)) for a, v in acc.items() if len(v) >= min_dom and a in PACE}
    if len(ctr) < 8: return np.nan, 0
    k = sorted(ctr)
    return float(stats.spearmanr([PACE[a] for a in k], [ctr[a] for a in k])[0]), len(k)

rep = {"note": "complete 748-domain set, same loader as the manuscript"}
for kind in ("natural", "designed"):
    sub = d[d["kind"] == kind]
    cell = cells(sub); doms = sorted(cell)
    r, n = rho_from(cell, doms)
    boot = []
    for _ in range(2000):
        rr, _ = rho_from(cell, list(RNG.choice(doms, len(doms), replace=True)))
        if np.isfinite(rr): boot.append(rr)
    rep[kind] = {"domains": int(sub["domain"].nunique()), "rows": int(len(sub)),
                 "n_residues": n, "spearman": round(r, 4),
                 "ci95": [round(float(q), 4) for q in np.quantile(boot, [.025, .975])]}
json.dump(rep, open(out, "w"), indent=2); print(json.dumps(rep, indent=2))
