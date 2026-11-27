"""Figures for PREPRINT_indel.md.

Fig 1  deletion cost against insertion cost, binned, natural and designed
Fig 2  asymmetry by secondary structure, paired within protein, both datasets
Fig 3  the asymmetry as a function of position tolerance
Fig 4  the failed strand-versus-helix claim, natural against designed

Usage: indel_env python 36_figures.py <lib3_csv> <results_dir> <afdir> <outdir>
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
exec(open(os.path.join(HERE, "34_scale.py")).read().split("def main")[0])  # noqa: S102

C = {"E": "#2f6f9f", "H": "#c2553a", "L": "#8a8f98"}
NAME = {"E": "beta strand", "H": "alpha helix", "L": "loop"}
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                     "figure.dpi": 160})


def pair_class(df, a, b, col="asym", nmin=2):
    return [float(g[g["ss"] == a][col].mean() - g[g["ss"] == b][col].mean())
            for _, g in df.groupby("domain")
            if (g["ss"] == a).sum() >= nmin and (g["ss"] == b).sum() >= nmin]


def ci(v, n=20000):
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    if len(v) < 5:
        return np.nan, np.nan, np.nan
    rng = np.random.default_rng(0)
    bs = v[rng.integers(0, len(v), size=(n, len(v)))].mean(axis=1)
    return v.mean(), np.quantile(bs, .025), np.quantile(bs, .975)


def main(lib3, rdir, afdir, outdir):
    os.makedirs(outdir, exist_ok=True)
    D, wtD = designed(lib3, afdir)
    N, wtN = natural(rdir)

    # ---- Fig 1
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.1))
    for ax, (nm, df, lim) in zip(axes, [("387 natural domains", N, 6.5),
                                        ("57 designed miniproteins", D, 2.2)]):
        d = df[df["ins"].notna() & df["del"].notna()].copy()
        d["b"] = pd.qcut(d["ins"], 12, duplicates="drop")
        gb = d.groupby("b", observed=True).agg(x=("ins", "mean"), y=("del", "mean"),
                                               n=("del", "size"))
        ax.plot([-1, lim], [-1, lim], color="#bbb", lw=1, ls="--", zorder=1)
        ax.plot(gb["x"], gb["y"], "o-", color="#2f6f9f", ms=4, lw=1.4, zorder=3)
        ax.set_xlabel("insertion cost (kcal/mol)")
        ax.set_ylabel("deletion cost (kcal/mol)")
        ax.set_title(nm, fontsize=9)
        ax.set_xlim(-0.6 if lim < 3 else -0.4, lim)
        ax.set_ylim(-0.6 if lim < 3 else -0.4, lim)
        ax.text(0.04, 0.93, "y = x", transform=ax.transAxes, color="#999", fontsize=8)
    fig.suptitle("Deletions cost more than insertions only where the position is tolerant",
                 fontsize=10, y=1.02)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "fig1_del_vs_ins.png"), bbox_inches="tight")
    plt.close(fig)

    # ---- Fig 2
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.0), sharey=False)
    for ax, (nm, df) in zip(axes, [("387 natural domains", N),
                                   ("57 designed miniproteins", D)]):
        for i, s in enumerate(["E", "H", "L"]):
            v = [float(g[g["ss"] == s]["asym"].mean()) for _, g in df.groupby("domain")
                 if (g["ss"] == s).sum() >= 2]
            m, lo, hi = ci(v)
            ax.errorbar(i, m, yerr=[[m - lo], [hi - m]], fmt="o", color=C[s], ms=6, capsize=4)
        ax.axhline(0, color="#ccc", lw=1)
        ax.set_xticks(range(3))
        ax.set_xticklabels([NAME[s] for s in ["E", "H", "L"]])
        ax.set_xlim(-0.5, 2.5)
        ax.set_ylabel("insertion cost - deletion cost\n(kcal/mol)")
        ax.set_title(nm, fontsize=9)
    fig.suptitle("Insertions are cheaper than deletions, most of all in beta strands",
                 fontsize=10, y=1.03)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "fig2_asymmetry_by_ss.png"), bbox_inches="tight")
    plt.close(fig)

    # ---- Fig 3
    fig, ax = plt.subplots(figsize=(4.2, 3.1))
    for nm, df, col in (("natural", N, "#2f6f9f"), ("designed", D, "#c2553a")):
        d = df[df["ins"].notna() & df["del"].notna()].copy()
        d["b"] = pd.qcut(d["ins"], 8, duplicates="drop")
        xs, ys, los, his = [], [], [], []
        for b, g in d.groupby("b", observed=True):
            v = [float((x["del"] - x["ins"]).mean()) for _, x in g.groupby("domain")
                 if len(x) >= 3]
            m, lo, hi = ci(v)
            if np.isfinite(m):
                xs.append(float(g["ins"].mean()))
                ys.append(m)
                los.append(lo)
                his.append(hi)
        ax.plot(xs, ys, "o-", color=col, ms=4, lw=1.4, label=nm)
        ax.fill_between(xs, los, his, color=col, alpha=0.15, lw=0)
    ax.axhline(0, color="#ccc", lw=1)
    ax.set_xlabel("insertion cost at that position (kcal/mol)")
    ax.set_ylabel("deletion cost - insertion cost\n(kcal/mol)")
    ax.set_title("The asymmetry needs headroom", fontsize=9)
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "fig3_headroom.png"), bbox_inches="tight")
    plt.close(fig)

    # ---- Fig 4
    fig, ax = plt.subplots(figsize=(4.4, 3.0))
    rows = [("insertion\nnatural", pair_class(N, "E", "H", "ins"), "#2f6f9f"),
            ("insertion\ndesigned", pair_class(D, "E", "H", "ins"), "#c2553a"),
            ("deletion\nnatural", pair_class(N, "E", "H", "del"), "#2f6f9f"),
            ("deletion\ndesigned", pair_class(D, "E", "H", "del"), "#c2553a")]
    for i, (lab, v, col) in enumerate(rows):
        m, lo, hi = ci(v)
        ax.errorbar(m, i, xerr=[[m - lo], [hi - m]], fmt="o", color=col, ms=6, capsize=4)
    ax.axvline(0, color="#ccc", lw=1)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([r[0] for r in rows], fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("strand minus helix cost (kcal/mol)")
    ax.set_title("A within-protein effect that did not replicate", fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "fig4_strand_vs_helix.png"), bbox_inches="tight")
    plt.close(fig)
    print("wrote 4 figures to", outdir)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
