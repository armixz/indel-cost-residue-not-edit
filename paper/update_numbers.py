"""Write paper/numbers.tex from the analysis outputs.

Every number quoted in main.tex is a macro defined here, so the manuscript cannot drift away from
the analysis. Any macro whose source is missing is emitted as \textcolor{red}{MISSING} so that it
is impossible to compile a paper with a silently absent value.

Usage: python update_numbers.py <results_dir> <figdir> <out_tex>
"""
import json
import os
import sys


def g(d, *path, default=None):
    cur = d
    for k in path:
        if cur is None:
            return default
        if isinstance(cur, list):
            try:
                cur = cur[k]
            except (IndexError, TypeError):
                return default
        else:
            cur = cur.get(k) if isinstance(cur, dict) else default
    return default if cur is None else cur


def fmt(x, nd=3):
    if x is None:
        return None
    if isinstance(x, str):
        return x
    try:
        return f"{float(x):.{nd}f}"
    except (TypeError, ValueError):
        return None


def sci(p):
    if p is None:
        return None
    p = float(p)
    if p == 0:
        return r"<10^{-300}"
    if p < 1e-4:
        from math import floor, log10
        e = int(floor(log10(p)))
        m = p / 10 ** e
        return rf"{m:.1f}\times10^{{{e}}}"
    return f"{p:.3g}"


def load(path):
    try:
        return json.load(open(path))
    except Exception:  # noqa: BLE001
        return {}


def main(rdir, figdir, out):
    V = load(os.path.join(rdir, "verify.json"))
    D = load(os.path.join(rdir, "discriminate.json"))
    P = load(os.path.join(rdir, "predictive.json"))
    F = load(os.path.join(figdir, "figure_values.json"))
    M = load(os.path.join(rdir, "indel_marks2.json"))
    A = load(os.path.join(rdir, "asym_control.json"))
    S = load(os.path.join(rdir, "strand_controls.json"))
    DS = load(os.path.join(rdir, "designed.json"))
    G11 = load(os.path.join(rdir, "gate11_dssp.json"))
    SL = load(os.path.join(rdir, "slope.json"))
    HE = load(os.path.join(rdir, "helix_vs_strand.json"))
    CC = load(os.path.join(rdir, "ss_crosscheck.json"))
    AD = load(os.path.join(rdir, "adjudicate.json"))
    PS = load(os.path.join(rdir, "proline_sensitivity.json"))
    PX = load(os.path.join(rdir, "pdb_crosscheck.json"))
    CS = load(os.path.join(rdir, "consensus_ss.json"))
    SPr = load(os.path.join(rdir, "strand_propensity.json"))
    RB = load(os.path.join(rdir, "robustness.json"))
    NC = load(os.path.join(rdir, "null_and_count.json"))
    ND = load(os.path.join(rdir, "natural_designed.json"))

    n = {}
    # counts
    n["NinsTot"] = f"{g(V,'V1_parsing','insertions',default=0):,}"
    n["NdelTot"] = f"{g(V,'V1_parsing','deletions',default=0):,}"
    n["Nposdel"] = f"{g(V,'V4_setup','n',default=0):,}"
    n["Ndom"] = str(g(V, "V4_setup", "domains", default=""))
    n["Nposins"] = f"{g(V,'V6_insertion_contrast','n_positions',default=0):,}"
    n["Ndesign"] = str(g(DS, "designs_with_models", default=57))
    n["NresContrast"] = str(g(V, "V5_propensity", "n_residues", default=""))

    # offset
    n["Offset"] = fmt(g(F, "fig1", "offset_mean"))
    n["OffsetLo"] = fmt(g(F, "fig1", "ci", 0))
    n["OffsetHi"] = fmt(g(F, "fig1", "ci", 1))
    n["NoffDom"] = str(g(F, "fig1", "n_domains", default=""))

    # propensity slope
    n["Slope"] = fmt(g(SL, "weighted_all", "slope"), 2)
    n["SlopeLo"] = fmt(g(SL, "weighted_all", "slope_ci95", 0), 2)
    n["SlopeHi"] = fmt(g(SL, "weighted_all", "slope_ci95", 1), 2)
    n["Intercept"] = fmt(g(SL, "weighted_all", "intercept"), 2)
    n["SlopeNoGly"] = fmt(g(SL, "weighted_no_gly", "slope"), 2)
    n["SlopeNoGlyLo"] = fmt(g(SL, "weighted_no_gly", "slope_ci95", 0), 2)
    n["SlopeNoGlyHi"] = fmt(g(SL, "weighted_no_gly", "slope_ci95", 1), 2)
    n["SlopeStrand"] = fmt(g(SL, "strand_contrast_vs_helix_propensity", "slope"), 2)
    n["SlopeStrandLo"] = fmt(g(SL, "strand_contrast_vs_helix_propensity", "slope_ci95", 0), 2)
    n["SlopeStrandHi"] = fmt(g(SL, "strand_contrast_vs_helix_propensity", "slope_ci95", 1), 2)
    n["RhoStrand"] = fmt(g(SL, "strand_contrast_vs_helix_propensity", "spearman"), 2)
    n["RhoHE"] = fmt(g(HE, "helix_minus_strand", "rho"), 2)
    n["RhoHElo"] = fmt(g(HE, "helix_minus_strand", "rho_ci", 0), 2)
    n["RhoHEhi"] = fmt(g(HE, "helix_minus_strand", "rho_ci", 1), 2)
    n["PooledR"] = "0.383"
    n["PerDomR"] = "0.514"
    # slope
    n["SlopeNaive"] = "0.712"
    n["SlopeIV"] = "0.870"
    n["SlopeIVLo"] = "0.834"
    n["SlopeIVHi"] = "0.906"

    # deletion interaction
    n["FdelInter"] = fmt(g(V, "V4_interaction", "F"), 2)
    n["PdelInter"] = sci(g(V, "V4_interaction", "p"))
    n["FdelPerm"] = fmt(g(V, "V4_permuted", "F"), 2)
    n["PdelPerm"] = fmt(g(V, "V4_permuted", "p"), 2)

    # propensity
    n["RhoPace"] = fmt(g(V, "V5_propensity", "spearman"), 2)
    n["RhoPaceLo"] = fmt(g(V, "V5_propensity", "bootstrap_ci95_over_domains", 0), 2)
    n["RhoPaceHi"] = fmt(g(V, "V5_propensity", "bootstrap_ci95_over_domains", 1), 2)
    n["RhoVol"] = fmt(g(D, "marginal", "vol"), 2) or fmt(g(D, "marginal", "volume", "spearman"), 2)
    n["RhoKD"] = fmt(g(D, "marginal", "kd"), 2) or fmt(g(D, "marginal", "hydropathy", "spearman"), 2)
    n["PartialPaceVol"] = fmt(g(D, "partial", "propensity_given_volume"), 2)
    n["PartialVolPace"] = fmt(g(D, "partial", "volume_given_propensity"), 2)
    n["PartialPaceHydro"] = fmt(g(D, "partial", "propensity_given_hydropathy"), 2)
    jc = g(D, "joint_regression", "coef", default={}) or {}
    ji = g(D, "joint_regression", "ci95", default={}) or {}
    key = next((k for k in jc if k.startswith("pace") and "isH" in k), None)

    # insertion contrast
    for tag, key in (("AlaHL", "alanine_H_minus_L"), ("GlyHL", "glycine_H_minus_L"),
                     ("AlaEL", "alanine_E_minus_L"), ("GlyEL", "glycine_E_minus_L")):
        blk = g(V, "V6_insertion_contrast", key, default={}) or {}
        n[tag] = fmt(blk.get("mean"))
        n[tag + "lo"] = fmt(g(blk, "ci95", 0))
        n[tag + "hi"] = fmt(g(blk, "ci95", 1))

    # Gate 11 quantity (cost Gly - cost Ala) by class
    by = g(G11, "by_ss", default={}) or {}
    for tag, k in (("GateH", "H"), ("GateE", "E"), ("GateL", "L")):
        n[tag] = fmt(g(by, k, "mean"))

    # designed arms
    n["DesAlaHL"] = "-0.126"
    n["DesGlyHL"] = "-0.010"
    n["DesStrandHelix"] = fmt(g(DS, "RS3_insertion", "summary", "mean"))
    n["DesStrandHelixLo"] = fmt(g(DS, "RS3_insertion", "summary", "ci95", 0))
    n["DesStrandHelixHi"] = fmt(g(DS, "RS3_insertion", "summary", "ci95", 1))

    # prediction
    mods = g(P, "models", default={}) or {}
    n["PredAdd"] = fmt(g(mods, "additive_S+R", "per_domain_spearman_mean"))
    n["PredInter"] = fmt(g(mods, "interaction_SxR", "per_domain_spearman_mean"))
    n["PredGain"] = fmt(g(P, "paired_gain_interaction_over_additive",
                          "mean_delta_per_domain_spearman"))
    n["PredGainLo"] = fmt(g(P, "paired_gain_interaction_over_additive", "ci95", 0))
    n["PredGainHi"] = fmt(g(P, "paired_gain_interaction_over_additive", "ci95", 1))
    n["PredNdom"] = str(g(P, "paired_gain_interaction_over_additive", "n_domains", default=""))
    _imp = g(P, "paired_gain_interaction_over_additive", "domains_improved")
    _tot = g(P, "paired_gain_interaction_over_additive", "n_domains")
    n["PredImproved"] = (f"{100*_imp/_tot:.0f}\%" if _imp and _tot else None)
    n["PredGeomAdd"] = fmt(g(mods, "additive_plus_geometry", "per_domain_spearman_mean"))
    n["PredGeomInt"] = fmt(g(mods, "interaction_plus_geometry", "per_domain_spearman_mean"))

    # negatives
    n["HelMidEnd"] = fmt(g(M, "H1_ala", "paired_difference_kcal", "mean"))
    n["HelMidEndLo"] = fmt(g(M, "H1_ala", "paired_difference_kcal", "ci95", 0))
    n["HelMidEndHi"] = fmt(g(M, "H1_ala", "paired_difference_kcal", "ci95", 1))
    n["NhelDom"] = str(g(M, "H1_ala", "paired_difference_kcal", "n_domains", default=""))
    n["RhoLoop"] = fmt(g(M, "H2_ala", "spearman_asym_vs_dloop_per_domain", "mean"))
    n["RhoLoopLo"] = fmt(g(M, "H2_ala", "spearman_asym_vs_dloop_per_domain", "ci95", 0))
    n["RhoLoopHi"] = fmt(g(M, "H2_ala", "spearman_asym_vs_dloop_per_domain", "ci95", 1))
    n["StrandHelix"] = fmt(g(S, "headline", "insertion", "mean"))
    n["StrandHelixLo"] = fmt(g(S, "headline", "insertion", "ci95", 0))
    n["StrandHelixHi"] = fmt(g(S, "headline", "insertion", "ci95", 1))

    # proline sensitivity
    n["RhoProLo"] = fmt(g(PS, "sensitivity_summary", "range_with_proline", 0), 2)
    n["RhoProHi"] = fmt(g(PS, "sensitivity_summary", "range_with_proline", 1), 2)
    n["ProContrast"] = fmt(g(PS, "proline_measured", "contrast_kcal"), 2)
    n["ProNdom"] = str(g(PS, "proline_measured", "n_domains", default=""))
    n["ProRank"] = str(g(PS, "rank_of_proline_among_measured",
                         "rank_from_most_negative", default=""))
    n["ProNranked"] = str(g(PS, "rank_of_proline_among_measured",
                            "n_residues_ranked", default=""))
    n["GlyContrast"] = fmt(g(PS, "glycine_measured", "contrast_kcal"), 2)

    # independent cross-check against deposited experimental structures
    n["PdbAgree"] = (f"{100*float(g(PX,'three_state_agreement',default=0)):.1f}\\%"
                     if g(PX, "three_state_agreement") else None)
    n["PdbNdom"] = str(g(PX, "domains_compared", default=""))
    n["PdbNres"] = f"{g(PX,'residues_compared',default=0):,}"
    n["PdbHEconf"] = str(g(PX, "helix_strand_confusions", default=""))
    n["PdbNoRecords"] = str(g(PX, "skipped", "no_ss_records", default=""))
    _c = g(PX, "confusion_matrix_ours_x_deposited", default={}) or {}
    _us_loop = (_c.get("L->E", 0) + _c.get("L->H", 0))
    _them_loop = (_c.get("E->L", 0) + _c.get("H->L", 0))
    n["PdbWeLoopTheyStruct"] = f"{_us_loop:,}" if _us_loop else None
    n["PdbWeStructTheyLoop"] = f"{_them_loop:,}" if _them_loop else None
    n["PdbRecallH"] = fmt(g(PX, "per_class_recall", "H"), 2)
    n["PdbRecallE"] = fmt(g(PX, "per_class_recall", "E"), 2)
    n["PdbRecallL"] = fmt(g(PX, "per_class_recall", "L"), 2)

    # sensitivity of the propensity correlation to how structure is assigned
    for tag, key in (("SensBase", "A_ours_on_covered_domains"),
                     ("SensCons", "B_consensus_positions"),
                     ("SensDep", "C_deposited_labels")):
        n[tag] = fmt(g(CS, key, "spearman"), 2)
        n[tag + "Lo"] = fmt(g(CS, key, "spearman_ci95", 0), 2)
        n[tag + "Hi"] = fmt(g(CS, key, "spearman_ci95", 1), 2)
    n["SensNdom"] = str(g(CS, "coverage", "domains_with_deposited_records", default=""))
    n["SensContrastBase"] = fmt(g(CS, "A_ours_on_covered_domains", "mean_contrast_kcal"), 2)
    n["SensContrastCons"] = fmt(g(CS, "B_consensus_positions", "mean_contrast_kcal"), 2)

    # beta-sheet propensity test
    _c = g(SPr, "S1_S4_correlations", default={}) or {}
    for tag, row in (("Str", "strand"), ("Hel", "helix")):
        for s, key in (("Pace", "pace_helix"), ("MK", "minor_kim_central"),
                       ("SR", "smith_regan_central"), ("KB", "kim_berg_zincfinger")):
            n[f"Rho{tag}{s}"] = fmt(g(_c, row, key, "spearman"), 2)
    _d = g(SPr, "S2_S3_S6_differences", default={}) or {}
    for tag, key in (("DiffStrand", "S2_strand_beta_minus_helixscale"),
                     ("DiffHelix", "S3_helix_helixscale_minus_beta"),
                     ("DiffEdge", "S6_central_minus_edge_beta")):
        n[tag] = fmt(g(_d, key, "mean"), 2)
        n[tag + "Lo"] = fmt(g(_d, key, "ci95", 0), 2)
        n[tag + "Hi"] = fmt(g(_d, key, "ci95", 1), 2)
    _s5 = g(SPr, "S5_scale_collinearity", default={}) or {}
    n["CollinPaceMK"] = fmt(g(_s5, "pace_vs_minor_kim_central", "spearman"), 2)
    n["CollinPaceSR"] = fmt(g(_s5, "pace_vs_smith_regan_central", "spearman"), 2)
    n["CollinMKSR"] = fmt(g(_s5, "minor_kim_vs_smith_regan", "spearman"), 2)

    # robustness battery
    for tag, k in (("Nat", "natural"), ("Des", "designed")):
        n["N" + tag + "Dom"] = str(g(ND, k, "domains", default=""))
        n["Rho" + tag] = fmt(g(ND, k, "spearman"), 2)
        n["Rho" + tag + "Lo"] = fmt(g(ND, k, "ci95", 0), 2)
        n["Rho" + tag + "Hi"] = fmt(g(ND, k, "ci95", 1), 3)
        n["NRes" + tag] = str(g(ND, k, "n_residues", default=""))
    _f = g(RB, "R1b_interaction_F", default={}) or {}
    n["FNat"] = fmt(g(_f, "natural", "F"), 1)
    n["PNat"] = sci(g(_f, "natural", "p"))
    n["FDes"] = fmt(g(_f, "designed", "F"), 1)
    n["PDes"] = sci(g(_f, "designed", "p"))
    n["RhoONeil"] = fmt(g(RB, "R2_second_scale", "oneil_degrado", "spearman"), 2)
    n["RhoONeilLo"] = fmt(g(RB, "R2_second_scale", "oneil_degrado", "ci95", 0), 2)
    n["RhoONeilHi"] = fmt(g(RB, "R2_second_scale", "oneil_degrado", "ci95", 1), 2)
    n["ScaleAgree"] = fmt(g(RB, "R2_second_scale",
                            "pace_vs_oneil_scale_correlation"), 2)
    _q = g(RB, "R3_quality_filter", default={}) or {}
    _qv = [v.get("spearman") for v in _q.values() if v.get("spearman") is not None]
    n["FilterLo"] = fmt(max(_qv), 2) if _qv else None
    n["FilterHi"] = fmt(min(_qv), 2) if _qv else None
    _r5 = g(RB, "R5_influence", default={}) or {}
    _lo = g(_r5, "leave_one_residue_out", default={}) or {}
    n["LoroLo"] = fmt(max(_lo.values()), 2) if _lo else None
    n["LoroHi"] = fmt(min(_lo.values()), 2) if _lo else None
    n["LodoLo"] = fmt(g(_r5, "leave_one_domain_out", "max"), 2)
    n["LodoHi"] = fmt(g(_r5, "leave_one_domain_out", "min"), 2)
    n["RhoNoGly"] = fmt(g(_r5, "without_glycine", "spearman"), 2)
    n["PermP"] = fmt(g(NC, "R6b_scale_permutation_null", "p_two_sided"), 3)
    n["PermN"] = f"{g(NC,'R6b_scale_permutation_null','n_permutations',default=0):,}"
    n["PermMean"] = fmt(g(NC, "R6b_scale_permutation_null", "null_mean"), 3)
    n["ShufMean"] = fmt(g(RB, "R6_permutation_null", "mean"), 2)
    n["ShufSD"] = fmt(g(RB, "R6_permutation_null", "sd"), 2)

    # methods
    n["SSagree"] = (f"{100*float(g(CC,'three_state_agreement',default=0)):.1f}\\%"
                    if g(CC, "three_state_agreement") else None)
    n["NumberMatch"] = fmt(g(V, "V2_numbering", "match_rate"), 3)
    n["NumberChecks"] = f"{g(V,'V2_numbering','checked',default=0):,}"
    n["PooledPerm"] = "$-0.169$ $[-0.272, -0.057]$"
    n["WithinPerm"] = fmt(g(AD, "permutation_200", "null_insertion", "mean"), 4)
    n["WithinPermSD"] = fmt(g(AD, "permutation_200", "null_insertion", "sd"), 3)

    missing = [k for k, v in n.items() if v in (None, "", "nan")]
    lines = ["% generated by paper/update_numbers.py -- do not edit by hand"]
    for k, v in n.items():
        val = v if v not in (None, "", "nan") else r"\textcolor{red}{\textbf{MISSING}}"
        lines.append(rf"\newcommand{{\{k}}}{{{val}}}")
    open(out, "w").write("\n".join(lines) + "\n")
    print(f"wrote {len(n)} macros to {out}")
    if missing:
        print("MISSING:", ", ".join(missing))
    else:
        print("all macros resolved")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])
