"""Publication figures for the indel paper. Vector PDF output for LaTeX.

Fig 1  the experiment: deletion vs insertion cost, and the mean offset
Fig 2  deletions, all 20 residues: helix-minus-loop contrast vs helix propensity
Fig 3  insertions, two residues: the structure contrast reverses between Ala and Gly
Fig 4  predictive gain from modelling the interaction

Usage: indel_env python make_figures.py <lib3_csv> <afdir> <results_dir> <outdir>
"""
import csv
import importlib.util
import json
import os
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ANA = os.path.join(os.path.dirname(HERE), "analysis")
sp = importlib.util.spec_from_file_location("d17", os.path.join(ANA, "17_dssp.py"))
D17 = importlib.util.module_from_spec(sp)
sp.loader.exec_module(D17)
PAT = re.compile(r"^(?P<dom>.+?)_(?P<kind>ins|del)(?P<res>[A-Z])(?P<pos>\d+)$")
RNG = np.random.default_rng(0)

# Pace & Scholtz 1998, Biophys J 75:422, abstract Table values, kcal/mol relative to Ala.
# Proline is deliberately absent: the published abstract excludes it and the circulating
# value of 3.16 could not be verified against the paper.
PACE = {"A": 0.00, "L": 0.21, "R": 0.21, "M": 0.24, "K": 0.26, "Q": 0.39, "E": 0.40,
        "I": 0.41, "W": 0.49, "S": 0.50, "Y": 0.53, "F": 0.54, "H": 0.61, "V": 0.61,
        "N": 0.65, "T": 0.66, "C": 0.68, "D": 0.69, "G": 1.00}

BLUE, RED, GREY = "#2f6f9f", "#c2553a", "#8a8f98"
plt.rcParams.update({
    "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 8.5, "legend.fontsize": 7,
    "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "axes.spines.top": False,
    "axes.spines.right": False, "axes.linewidth": 0.7, "xtick.major.width": 0.7,
    "ytick.major.width": 0.7, "pdf.fonttype": 42, "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans"]})


def load(lib3, afdir, rdir):
    csv.field_size_limit(10 ** 7)
    wt, rows = {}, []
    with open(lib3, newline="") as fh:
        for r in csv.DictReader(fh):
            nm = r["name"]
            if nm.rsplit("_", 1)[-1].startswith("wt"):
                wt.setdefault(nm.rsplit("_", 1)[0], []).append(float(r["deltaG"]))
                continue
            m = PAT.match(nm)
            if not m:
                continue
            try:
                rows.append({"domain": m["dom"], "kind": m["kind"], "res": m["res"],
                             "pos": int(m["pos"]), "dG": float(r["deltaG"]),
                             "ci": float(r["deltaG_95CI"])})
            except (ValueError, KeyError):
                continue
    wtm = {k: float(np.mean(v)) for k, v in wt.items()}
    v = pd.DataFrame(rows)
    v = v[v["ci"] <= 1.0].copy()
    v["wtdG"] = v["domain"].map(wtm)
    v = v[v["wtdG"].notna()]
    v["cost"] = v["wtdG"] - v["dG"]
    ss = []
    for dom in sorted(v["domain"].unique()):
        p = os.path.join(afdir, dom.replace(".pdb", "") + ".pdb")
        if not os.path.exists(p):
            continue
        res = D17.parse_backbone(p)
        s, keys = D17.assign(res)
        if not s:
            continue
        el = D17.elements(s, keys)
        for i in keys:
            ss.append({"domain": dom, "pos": i, "aa": res[i]["aa"], "ss": el[i][0]})
    S = pd.DataFrame(ss)
    for key, kind, res in (("insA", "ins", "A"), ("insG", "ins", "G"), ("del", "del", None)):
        sub = v[v["kind"] == kind]
        if res:
            sub = sub[sub["res"] == res]
        S = S.merge(sub.groupby(["domain", "pos"])["cost"].mean().rename(key).reset_index(),
                    on=["domain", "pos"], how="left")
    return S, wtm


def ci(v, n=20000):
    v = np.asarray([x for x in v if np.isfinite(x)], float)
    if len(v) < 5:
        return np.nan, np.nan, np.nan
    bs = v[RNG.integers(0, len(v), size=(n, len(v)))].mean(axis=1)
    return v.mean(), np.quantile(bs, .025), np.quantile(bs, .975)


def contrast(df, col, a, b, nmin=2):
    return [float((x[x["ss"] == a][col] - x[x["ss"] == a]["del"]).mean()
                  - (x[x["ss"] == b][col] - x[x["ss"] == b]["del"]).mean())
            for _, x in df.groupby("domain")
            if (x["ss"] == a).sum() >= nmin and (x["ss"] == b).sum() >= nmin]


def fig1(S, out):
    d = S[S["insA"].notna() & S["del"].notna()].copy()
    fig, ax = plt.subplots(1, 2, figsize=(5.6, 2.4))
    d["b"] = pd.qcut(d["insA"], 14, duplicates="drop")
    gb = d.groupby("b", observed=True).agg(x=("insA", "mean"), y=("del", "mean"))
    ax[0].plot([-1, 6.5], [-1, 6.5], color="#ccc", lw=0.8, ls="--", zorder=1)
    ax[0].plot(gb["x"], gb["y"], "o-", color=BLUE, ms=3, lw=1.1, zorder=3)
    ax[0].set_xlim(-0.6, 6.5)
    ax[0].set_ylim(-0.6, 6.5)
    ax[0].set_xlabel(r"insertion cost $\Delta\Delta G$ (kcal mol$^{-1}$)")
    ax[0].set_ylabel(r"deletion cost (kcal mol$^{-1}$)")
    ax[0].set_title("a   deletions cost more", loc="left", fontweight="bold")
    ax[0].text(3.4, 2.4, "y = x", color="#999", fontsize=7)

    off = [float((g["del"] - g["insA"]).mean()) for _, g in d.groupby("domain") if len(g) >= 5]
    m, lo, hi = ci(off)
    ax[1].hist(off, bins=40, color=BLUE, alpha=.8, lw=0)
    ax[1].axvline(0, color="#999", lw=0.8)
    ax[1].axvline(m, color=RED, lw=1.2)
    ax[1].set_xlabel(r"per-domain mean offset (kcal mol$^{-1}$)")
    ax[1].set_ylabel("domains")
    ax[1].set_title("b   offset across domains", loc="left", fontweight="bold")
    ax[1].text(.97, .93, f"mean {m:.3f}\n95% CI [{lo:.3f}, {hi:.3f}]\n{len(off)} domains",
               transform=ax[1].transAxes, ha="right", va="top", fontsize=6.8, color="#444")
    fig.tight_layout()
    fig.savefig(os.path.join(out, "fig1.pdf"), bbox_inches="tight")
    plt.close(fig)
    return {"offset_mean": round(float(m), 4), "ci": [round(float(lo), 4), round(float(hi), 4)],
            "n_domains": len(off)}


def fig2(S, out):
    d = S[S["del"].notna()].copy()
    pts = []
    for aa, g in d.groupby("aa"):
        if aa not in PACE:
            continue
        v = [float(x[x["ss"] == "H"]["del"].mean() - x[x["ss"] == "L"]["del"].mean())
             for _, x in g.groupby("domain")
             if (x["ss"] == "H").sum() >= 1 and (x["ss"] == "L").sum() >= 1]
        if len(v) >= 5:
            m, lo, hi = ci(v)
            pts.append((aa, PACE[aa], m, lo, hi, len(v)))
    x = np.array([p[1] for p in pts])
    y = np.array([p[2] for p in pts])
    rho = stats.spearmanr(x, y)[0]
    fig, ax = plt.subplots(figsize=(3.3, 2.7))
    for aa, px, py, lo, hi, n in pts:
        ax.errorbar(px, py, yerr=[[py - lo], [hi - py]], fmt="none", ecolor="#ccd",
                    elinewidth=0.7, zorder=2)
        ax.scatter(px, py, s=10 + n / 3, color=BLUE, zorder=3, alpha=.9, edgecolor="none")
        ax.annotate(aa, (px, py), textcoords="offset points", xytext=(3.5, 2.5), fontsize=6.5,
                    color="#333")
    b, a = np.polyfit(x, y, 1)
    xs = np.linspace(x.min() - .03, x.max() + .03, 50)
    ax.plot(xs, a + b * xs, color=RED, lw=1.1, zorder=1)
    ax.axhline(0, color="#ddd", lw=0.7)
    ax.set_xlabel(r"helix propensity of the deleted residue" "\n"
                  r"(Pace \& Scholtz, kcal mol$^{-1}$ vs Ala)")
    ax.set_ylabel("deletion cost, helix $-$ loop\n" r"(kcal mol$^{-1}$)")
    ax.text(.96, .95, rf"Spearman $\rho$ = {rho:.2f}", transform=ax.transAxes, ha="right",
            va="top", fontsize=7, color="#444")
    fig.tight_layout()
    fig.savefig(os.path.join(out, "fig2.pdf"), bbox_inches="tight")
    plt.close(fig)
    return {"spearman": round(float(rho), 4), "n_residues": len(pts),
            "contrasts": {p[0]: round(float(p[2]), 4) for p in pts}}


def fig3(S, out):
    t = S[S["insA"].notna() & S["insG"].notna() & S["del"].notna()].copy()
    rows, res = [], {}
    for col, lab, c in (("insA", "alanine", BLUE), ("insG", "glycine", RED)):
        for a, b, nm in (("H", "L", "helix $-$ loop"), ("E", "L", "strand $-$ loop")):
            v = contrast(t, col, a, b)
            m, lo, hi = ci(v)
            rows.append((f"{nm}\n{lab} insertion", m, lo, hi, c))
            res[f"{lab}_{a}_{b}"] = [round(float(m), 4), round(float(lo), 4), round(float(hi), 4),
                                     len(v)]
    fig, ax = plt.subplots(figsize=(3.4, 2.5))
    yy = np.arange(len(rows))
    for i, (lab, m, lo, hi, c) in enumerate(rows):
        ax.errorbar(m, yy[i], xerr=[[m - lo], [hi - m]], fmt="o", color=c, ms=4.5, capsize=2.5,
                    elinewidth=1)
    ax.axvline(0, color="#bbb", lw=0.8)
    ax.set_yticks(yy)
    ax.set_yticklabels([r[0] for r in rows], fontsize=6.8)
    ax.invert_yaxis()
    ax.set_xlabel("asymmetry relative to loops\n" r"(insertion $-$ deletion, kcal mol$^{-1}$)")
    fig.tight_layout()
    fig.savefig(os.path.join(out, "fig3.pdf"), bbox_inches="tight")
    plt.close(fig)
    return res


def fig4(rdir, out):
    p = os.path.join(rdir, "predictive.json")
    if not os.path.exists(p):
        return {"status": "predictive.json missing"}
    r = json.load(open(p))
    order = ["structure_only", "residue_only", "additive_S+R", "interaction_SxR",
             "additive_plus_geometry", "interaction_plus_geometry"]
    names = ["structure", "residue", "structure + residue", "structure $\\times$ residue",
             "+ geometry", "+ geometry, interaction"]
    vals = [r["models"][k]["per_domain_spearman_mean"] for k in order if k in r["models"]]
    names = names[:len(vals)]
    fig, ax = plt.subplots(figsize=(3.4, 2.4))
    cols = [GREY, GREY, GREY, BLUE, GREY, BLUE][:len(vals)]
    ax.barh(np.arange(len(vals)), vals, color=cols, height=.62)
    ax.set_yticks(np.arange(len(vals)))
    ax.set_yticklabels(names, fontsize=6.8)
    ax.invert_yaxis()
    ax.set_xlabel("per-domain Spearman, leave-one-domain-out")
    for i, v in enumerate(vals):
        ax.text(v + .004, i, f"{v:.3f}", va="center", fontsize=6.5, color="#444")
    ax.set_xlim(0, max(vals) * 1.18)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "fig4.pdf"), bbox_inches="tight")
    plt.close(fig)
    return {k: r["models"][k]["per_domain_spearman_mean"] for k in order if k in r["models"]}


def figS1(S, out):
    """The deletion penalty against position tolerance, binned two ways.

    Binning on the same alanine insertion cost that enters the offset produces a spurious decline
    through zero, because measurement error on the binning axis is shared with the quantity
    plotted. The glycine insertion at the same position is an independent measurement of the same
    tolerance and does not have that defect.
    """
    m = S[S["insA"].notna() & S["insG"].notna() & S["del"].notna()].copy()
    fig, ax = plt.subplots(figsize=(3.6, 2.6))
    for binc, lab, col in (("insA", "binned on the same measurement", "#bbbbbb"),
                           ("insG", "binned on an independent measurement", BLUE)):
        d = m.copy()
        d["b"] = pd.qcut(d[binc], 8, duplicates="drop")
        xs, ys, lo, hi = [], [], [], []
        for _, g in d.groupby("b", observed=True):
            v = [float((x["del"] - x["insA"]).mean()) for _, x in g.groupby("domain")
                 if len(x) >= 3]
            mm, a, b = ci(v)
            if np.isfinite(mm):
                xs.append(float(g[binc].mean()))
                ys.append(mm)
                lo.append(a)
                hi.append(b)
        ax.plot(xs, ys, "o-", color=col, ms=3.2, lw=1.2, label=lab)
        ax.fill_between(xs, lo, hi, color=col, alpha=.15, lw=0)
    ax.axhline(0, color="#ccc", lw=0.8)
    ax.set_xlabel(r"insertion cost at that position (kcal mol$^{-1}$)")
    ax.set_ylabel("deletion $-$ insertion cost\n" r"(kcal mol$^{-1}$)")
    ax.legend(frameon=False, fontsize=6.5, loc="upper right")
    fig.tight_layout()
    fig.savefig(os.path.join(out, "figS1.pdf"), bbox_inches="tight")
    plt.close(fig)
    return {"note": "two binnings"}


def main(lib3, afdir, rdir, out):
    os.makedirs(out, exist_ok=True)
    S, wtm = load(lib3, afdir, rdir)
    rep = {"fig1": fig1(S, out), "fig2": fig2(S, out), "fig3": fig3(S, out),
           "fig4": fig4(rdir, out), "figS1": figS1(S, out)}
    json.dump(rep, open(os.path.join(out, "figure_values.json"), "w"), indent=2)
    print(json.dumps(rep, indent=2))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
