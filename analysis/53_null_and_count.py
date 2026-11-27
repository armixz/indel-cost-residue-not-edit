"""Two follow-ups the robustness battery raised, run before anything is written.

(1) The battery reports 691 domains where the manuscript reports 748. The two loaders
    differ: 43_discriminate recomputes secondary structure from the AlphaFold models,
    while 52_robustness merges the cached tsub_ss.csv. The discrepancy must be explained
    before either number is quoted.

(2) The residue-shuffling null in R6 is not centred at zero but at about +0.19. That is
    the opposite sign from the observed effect, so it cannot manufacture the finding, but
    it does mean the pre-registered wording ("the null must be centred near zero") is not
    satisfied as written. A cleaner null for the specific question "is this ORDERING
    special" is to permute the propensity values across residues, holding the measured
    contrasts fixed. That null is exact and has no estimator bias at all.
"""
import importlib.util, json, os, sys
import numpy as np, pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
def mod(name, fn):
    sp = importlib.util.spec_from_file_location(name, os.path.join(HERE, fn))
    m = importlib.util.module_from_spec(sp); sp.loader.exec_module(m); return m
D43 = mod("d43", "43_discriminate.py")
R52 = mod("r52", "52_robustness.py")

lib3, afdir, rdir, out = sys.argv[1:5]
rep = {}

a = D43.load(lib3, afdir, rdir)
b = R52.load_raw(lib3, afdir, rdir, 1.0)
A, B = set(a["domain"].unique()), set(b["domain"].unique())
rep["domain_count"] = {
    "dssp_recomputed_43": len(A), "cached_ss_table_52": len(B),
    "in_43_not_52": sorted(A - B)[:8], "n_in_43_not_52": len(A - B),
    "in_52_not_43": sorted(B - A)[:8], "n_in_52_not_43": len(B - A),
    "kinds_of_the_difference": {k: int(v) for k, v in
        pd.Series([R52.kind(x) for x in sorted(A - B)]).value_counts().items()}
        if A - B else {}}

ctr = R52.contrasts(b)
aas = sorted(x for x in ctr if x in R52.PACE)
y = np.array([ctr[x] for x in aas]); x = np.array([R52.PACE[x] for x in aas])
obs = float(stats.spearmanr(x, y)[0])
rng = np.random.default_rng(0)
null = [float(stats.spearmanr(rng.permutation(x), y)[0]) for _ in range(20000)]
rep["R6b_scale_permutation_null"] = {
    "description": "permute the propensity values across residues, contrasts held fixed",
    "n_permutations": len(null), "observed": round(obs, 4),
    "null_mean": round(float(np.mean(null)), 4), "null_sd": round(float(np.std(null)), 4),
    "p_two_sided": round(float(np.mean(np.abs(null) >= abs(obs))), 5),
    "n_residues": len(aas)}
json.dump(rep, open(out, "w"), indent=2)
print(json.dumps(rep, indent=2))
