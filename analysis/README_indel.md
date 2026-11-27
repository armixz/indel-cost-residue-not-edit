# Indel work, night of 25-26 Nov 2026 -- script index

Environment: `/root/envs/indel/bin/python` in WSL (pandas, scipy, statsmodels, pyreadr,
matplotlib). R for the Lehner lab's .rds files: `/root/micromamba/envs/r/bin/Rscript`.
Run from `analysis/`; paths in the examples are relative to it.

| script | what it does |
|---|---|
| 17_dssp.py | Kabsch-Sander secondary structure on the AlphaFold models -> results/tsub_ss.csv |
| 18_ss_check.py | validates 17 against dihedral geometry it never used |
| 19_gate11_dssp.py | recomputes Gate 11 with the corrected classes |
| 21_indel_marks2.py | marks B, H1, H2, H3 of PREREG_indel.md, within-domain estimators |
| 22_asym_control.py | residue identity, volume, burial controls on the asymmetry |
| 23_ss_crosscheck.py | agreement with the Lehner lab's own structure table |
| 25_bridges.py | beta-bridge partners, parallel/antiparallel, edge vs interior |
| 26_mechanism.py | marks M1-M3 of PREREG_bulge.md plus kill conditions |
| 27_within_domain.py | within-protein vs between-protein decomposition |
| 28_decompose.py | insertion and deletion costs relative to each domain's own loops |
| 29_multiplicative.py | rules out a pure scale effect |
| 30_strand_controls.py | C1-C5 of PREREG_strand.md and the aPCA replication |
| 31_adjudicate.py | 200-draw permutation, shared domains, dynamic range |
| 32_unfolded_check.py | raw K50 and the assay's unfolded-state model |
| 33_designed.py | the 57 designed HEEH_KT miniproteins as a replication set |
| 34_scale.py | natural and designed on one scale (also a module other scripts import) |
| 35_law.py | deletion-vs-insertion fit -- SUPERSEDED, see the correction in RESULTS_indel.md |
| 36_figures.py | first figure set (fig3 there is the biased version, kept for the record) |
| 37_decomp.py | side-chain decomposition, PREREG_decomposition.md (not run to completion) |
| 38_fig_corrected.py | corrected fig3 (independent binning) and fig5 (Ala vs Gly) |

Read `RESULTS_indel.md` first; it carries its own correction section for the two claims that
were withdrawn during the night.
