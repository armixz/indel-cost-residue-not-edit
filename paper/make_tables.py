"""LaTeX table fragments for the manuscript, generated from the analysis outputs.

tableS1.tex  per-residue helix-minus-loop and strand-minus-loop deletion contrasts, with the
             helix propensity used, the number of contributing proteins, and bootstrap intervals
tableS2.tex  model comparison for the prediction analysis

Usage: python make_tables.py <results_dir> <outdir>
"""
import json
import os
import sys

AA3 = {"A": "Ala", "C": "Cys", "D": "Asp", "E": "Glu", "F": "Phe", "G": "Gly", "H": "His",
       "I": "Ile", "K": "Lys", "L": "Leu", "M": "Met", "N": "Asn", "P": "Pro", "Q": "Gln",
       "R": "Arg", "S": "Ser", "T": "Thr", "V": "Val", "W": "Trp", "Y": "Tyr"}
PACE = {"A": 0.00, "L": 0.21, "R": 0.21, "M": 0.24, "K": 0.26, "Q": 0.39, "E": 0.40,
        "I": 0.41, "W": 0.49, "S": 0.50, "Y": 0.53, "F": 0.54, "H": 0.61, "V": 0.61,
        "N": 0.65, "T": 0.66, "C": 0.68, "D": 0.69, "G": 1.00}


def load(p):
    try:
        return json.load(open(p))
    except Exception:  # noqa: BLE001
        return {}


def main(rdir, out):
    os.makedirs(out, exist_ok=True)
    V = load(os.path.join(rdir, "verify.json"))
    W = load(os.path.join(rdir, "which_term.json"))
    P = load(os.path.join(rdir, "predictive.json"))

    # ---------------- Table S1
    con = (V.get("V5_propensity") or {}).get("contrasts") or {}
    hel = ((W.get("relative_to_domain_mean") or {}).get("helix") or {}).get("values") or {}
    loo = ((W.get("relative_to_domain_mean") or {}).get("loop") or {}).get("values") or {}
    rows = []
    for aa in sorted(con, key=lambda a: PACE.get(a, 99)):
        if aa not in PACE:
            continue
        val, ndom = con[aa][0], con[aa][1]
        rows.append((aa, PACE[aa], val, ndom, hel.get(aa), loo.get(aa)))

    lines = [
        r"\begin{table}[t]", r"\centering", r"\small",
        r"\caption{\textbf{Per-residue deletion contrasts.} For each amino acid, the mean cost of "
        r"deleting it from a helix minus the mean cost of deleting it from a loop, paired within "
        r"each protein and averaged over the proteins that contain both. The last two columns "
        r"split the contrast into its parts, each expressed relative to that protein's own mean "
        r"deletion cost, showing that the ordering by propensity sits in the helix term. "
        r"Residues are sorted by helix propensity. Proline is absent: the "
        r"\citet{pace1998} scale excludes it.}",
        r"\label{tab:perres}",
        r"\begin{tabular}{llrrrr}", r"\toprule",
        r"Residue & & Helix propensity & Helix $-$ loop & Helix term & Loop term \\",
        r" & & (kcal\,mol$^{-1}$) & (kcal\,mol$^{-1}$) & & \\",
        r"\midrule"]
    for aa, pace, val, ndom, h, l in rows:
        f = lambda x: "--" if x is None else f"{x:+.2f}"
        lines.append(rf"{AA3[aa]} & {aa} & {pace:.2f} & {val:+.2f} & {f(h)} & {f(l)} \\")
    lines += [r"\bottomrule", r"\end{tabular}", r"\end{table}"]
    open(os.path.join(out, "tableS1.tex"), "w").write("\n".join(lines) + "\n")

    # ---------------- Table S2
    mods = P.get("models") or {}
    order = [("structure_only", "secondary structure only"),
             ("residue_only", "residue identity only"),
             ("additive_S+R", "structure $+$ residue"),
             ("interaction_SxR", r"structure $\times$ residue"),
             ("additive_plus_geometry", r"structure $+$ residue $+$ geometry"),
             ("interaction_plus_geometry", r"structure $\times$ residue $+$ geometry"),
             ("additive_plus_geometry_x_structure",
              r"structure $+$ residue $+$ geometry $\times$ structure")]
    lines = [
        r"\begin{table}[t]", r"\centering", r"\small",
        r"\caption{\textbf{Leave-one-domain-out prediction of deletion cost.} Ridge models on the "
        r"same rows. Each domain is predicted by a model fitted without it. Per-domain Spearman "
        r"is the mean over held-out domains; the lower decile is reported because a mean can hide "
        r"a long tail.}",
        r"\label{tab:pred}",
        r"\begin{tabular}{lrrrr}", r"\toprule",
        r"Model & Features & RMSE & Per-domain $\rho$ & Lower decile \\", r"\midrule"]
    for key, label in order:
        m = mods.get(key)
        if not m:
            continue
        lines.append(rf"{label} & {m['n_features']} & {m['rmse']:.3f} & "
                     rf"{m['per_domain_spearman_mean']:.3f} & "
                     rf"{m['per_domain_lower_decile']:.3f} \\")
    lines += [r"\bottomrule", r"\end{tabular}", r"\end{table}"]
    open(os.path.join(out, "tableS2.tex"), "w").write("\n".join(lines) + "\n")
    print(f"wrote tableS1 ({len(rows)} residues) and tableS2 ({len(mods)} models)")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
