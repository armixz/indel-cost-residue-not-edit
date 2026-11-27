"""Corrected figures: independent binning, and the alanine/glycine comparison.

fig3 in 36_figures.py binned the offset (deletion minus insertion) by the insertion cost itself.
That is the same measurement on both axes, so measurement error alone produces a declining line
crossing zero near the mean. Binning instead by the GLYCINE insertion cost -- an independent
measurement of the same position -- removes that. fig5 shows that the beta-strand excess depends
on which residue is inserted.
"""
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                     "figure.dpi": 160})
RNG = np.random.default_rng(0)


def ci(v, n=20000):
    v = np.asarray([x for x in v if np.isfinite(x)], float)
    if len(v) < 5:
        return np.nan, np.nan, np.nan
    bs = v[RNG.integers(0, len(v), size=(n, len(v)))].mean(axis=1)
    return v.mean(), np.quantile(bs, .025), np.quantile(bs, .975)


def load(rdir):
    ind = pd.read_csv(f"{rdir}/tsuboyama_indels.csv")
    feat = pd.read_csv(f"{rdir}/tsub_features.csv")
    ss = pd.read_csv(f"{rdir}/tsub_ss.csv")
    ss["domain"] = ss["domain"].str.replace(".pdb", "", regex=False)
    wt = json.load(open(f"{rdir}/tsub_wt_dG.json"))
    ind = ind[ind["dG_ci"] <= 1.0].copy()
    ind["wt"] = ind["domain"].map(wt)
    ind = ind[ind["wt"].notna()]
    ind["cost"] = ind["wt"] - ind["dG"]
    df = feat[["domain", "pos", "struct_pos"]].merge(
        ss[["domain", "struct_pos", "ss"]], on=["domain", "struct_pos"])
    for key, kind, res in (("insA", "ins", "A"), ("insG", "ins", "G"), ("del", "del", None)):
        s = ind[ind["kind"] == kind]
        if res:
            s = s[s["res"] == res]
        df = df.merge(s.groupby(["domain", "pos"])["cost"].mean().rename(key).reset_index(),
                      on=["domain", "pos"], how="left")
    return df.dropna(subset=["insA", "insG", "del"]).copy()


def main(rdir, outdir):
    m = load(rdir)
    os.makedirs(outdir, exist_ok=True)

    # fig3 corrected
    fig, ax = plt.subplots(figsize=(4.6, 3.2))
    for binc, lab, col in (("insA", "binned by the same measurement (biased)", "#bbbbbb"),
                           ("insG", "binned by an independent measurement", "#2f6f9f")):
        d = m.copy()
        d["b"] = pd.qcut(d[binc], 8, duplicates="drop")
        xs, ys, lo, hi = [], [], [], []
        for b, g in d.groupby("b", observed=True):
            v = [float((x["del"] - x["insA"]).mean()) for _, x in g.groupby("domain")
                 if len(x) >= 3]
            mm, a, bb = ci(v)
            if np.isfinite(mm):
                xs.append(float(g[binc].mean()))
                ys.append(mm)
                lo.append(a)
                hi.append(bb)
        ax.plot(xs, ys, "o-", color=col, ms=4, lw=1.5, label=lab)
        ax.fill_between(xs, lo, hi, color=col, alpha=.15, lw=0)
    ax.axhline(0, color="#ccc", lw=1)
    ax.set_xlabel("insertion cost at that position (kcal/mol)")
    ax.set_ylabel("deletion cost - insertion cost\n(kcal/mol)")
    ax.set_title("The deletion penalty shrinks but does not vanish", fontsize=9)
    ax.legend(frameon=False, fontsize=7.5, loc="upper right")
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "fig3_headroom_corrected.png"), bbox_inches="tight")
    plt.close(fig)

    # fig5 alanine vs glycine
    fig, ax = plt.subplots(figsize=(4.6, 3.0))
    labels, means, los, his, cols = [], [], [], [], []
    for res, col in (("insA", "#2f6f9f"), ("insG", "#c2553a")):
        for a, b in (("E", "L"), ("H", "L")):
            v = [float((x[x["ss"] == a][res] - x[x["ss"] == a]["del"]).mean()
                       - (x[x["ss"] == b][res] - x[x["ss"] == b]["del"]).mean())
                 for _, x in m.groupby("domain")
                 if (x["ss"] == a).sum() >= 2 and (x["ss"] == b).sum() >= 2]
            mm, lo, hi = ci(v)
            labels.append(f"{'strand' if a=='E' else 'helix'} - loop\n{'Ala' if res=='insA' else 'Gly'} insertion")
            means.append(mm)
            los.append(lo)
            his.append(hi)
            cols.append(col)
    y = np.arange(len(means))
    for i in range(len(means)):
        ax.errorbar(means[i], y[i], xerr=[[means[i] - los[i]], [his[i] - means[i]]],
                    fmt="o", color=cols[i], ms=6, capsize=4)
    ax.axvline(0, color="#ccc", lw=1)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=7.5)
    ax.invert_yaxis()
    ax.set_xlabel("extra asymmetry relative to loops (kcal/mol)")
    ax.set_title("The strand excess depends on which residue is inserted", fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "fig5_ala_vs_gly.png"), bbox_inches="tight")
    plt.close(fig)
    print("wrote corrected figures")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
