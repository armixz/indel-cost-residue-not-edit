# Literature reconnaissance: indels, done before any analysis (26 Nov 2026)

Five parallel searches. Everything below is what the field already owns, so it is not available to
be rediscovered. Items marked *snippet* were not read in full text by the searcher.

## The experimental landscape

- **Tsuboyama et al. 2023** (Nature; Zenodo 7992926). Our data. Lib3 holds 95,184 single-residue
  insertions (Gly and Ala) and 46,660 single-residue deletions over ~400 domains, with AlphaFold
  models, in kcal/mol.
- **Topolska, Beltran & Lehner 2025**, Nat Commun 16:2617 (bioRxiv 2023.10.06.561180). The main
  re-analysis of these indels; filtered to 181 domains / 6,288 positions, 171 domains used. Adds
  new yeast abundancePCA on 9 domains with 1-3 aa indels, and two domains with all 20 amino acids
  inserted per position. Key published results:
  - deletions worse overall: 68.8% strongly reduce abundance vs 62.1% for insertions
  - loops more tolerant than structured regions (Cohen's d -0.194 del, -0.184 ins); strands slightly
    worse than helices
  - substitution scans show 3-4 periodicity in helices; **indels show no periodicity**
  - solvent accessibility predicts substitutions well (r 0.65-0.76) but indels poorly
  - substitution tolerance predicts indel tolerance in loops (r ~0.51-0.54) but not in structure
    (r 0.22 ins, 0.33 del)
  - 1-aa insertion vs 1-aa deletion correlation **r = 0.484** (their new yeast data, not Tsuboyama);
    r falls to 0.279 at 2 aa and 0.145 at 3 aa
  - **"both ends of alpha-helices are indeed slightly more tolerant"** (Fig 5d) — no kcal/mol given
  - **"proteins can locally accommodate insertions ... by pushing the indel out as a 'bulge' ...
    Our analyses suggest that this is a rare occurrence"**
  - INDELi predictor: per-domain **mode** r 0.620 (del) / 0.608 (ins); IQR 0.512-0.686 and
    0.393-0.646. Comparators: ESM1b 0.575 del, PROVEAN 0.428 ins / 0.324 del, CADD ~0.
  - They do **not** use pLDDT.
- **Gutierrez & Rocklin**, bioRxiv 2024.12.18.629072. The 103 stabilizing indels (ddG > 1) of the
  same data, in 48 wild types. Deleted sites are more backbone-strained (+0.33 REU; strand +0.50).
  Gly/Ala comparison only at the 36 sites an insertion stabilises. **This is the nearest prior work
  to our Gate 11** and it is confined to the stabilizing subset.
- **Vetter et al. 1996**, Protein Sci 5:2399 and **Heinz et al. 1994**, J Mol Biol 236:869. T4
  lysozyme: 32 insertions, 5 deletions, nine helices, 8 structures. Class averages 2.5 / 2.6 / 4.9
  kcal/mol (translocation / end looping-out / body bulge); deletions 3-5 kcal/mol. Mechanism:
  register shifts toward the N terminus, spilling residues into the preceding loop. Verified from
  the abstracts verbatim; full texts could not be opened. In extreme cases the shift propagates
  3-6 A over 20+ residues and even disrupts the hydrophobic core to preserve the helix.
- **Arpino et al. 2014**, PLoS One (PMC4058518). EGFP G4 deletion: K3 rotates ~120 deg, register
  restored from L7 onward; variant destabilized by 5.19 kcal/mol. One deletion, a 3-10 helix at a
  terminus. No energy attributed to the shift itself.
- Other indel DMS: Kir2.1 DIMPLE (Genome Biol 2023), Enterovirus A71 proteome scan
  (Nat Microbiol 2024; ~45,000 insertions, ~6,000 deletions; indels essentially excluded from
  helices, attributed to register change *as an assertion, not a test*), yeast DHFR folding sensor.

## Evolutionary side

- **Yang et al. 2024** (PMC11102076). 11,444 alignments, five mammalian genomes, AF2 secondary
  structure. Observed/expected indel ratios: strand 0.43, turn 0.55, helix 0.57, bend 0.59,
  no structure 1.55. **Insertions more depleted than deletions in 24 of 28 comparisons**; strand
  insertion O/E 0.26-0.39 (rodent) vs 0.41-0.52 (primate).
- MBE 2011, 28:291: structural displacement around an indel decays exponentially with distance, in
  homolog structure pairs. Geometry, not fitness.
- Beta bulges: Chan, Mazzarella & Thornton 1993, Protein Sci 2:1574 (362 bulges in 170 proteins);
  Acta Cryst D 2021 (classic 5133, G1alpha 3348, G1beta 1506). **No energetic cost of a bulge found
  anywhere.**
- Pi-helices as absorbed single-residue insertions: Cooley, Arp & Karplus 2010, JMB 404:232 — ~15%
  of proteins contain one, >95% unannotated, most arose from alpha-helices by a single insertion
  (*comparative evidence*). Counter: PMC6295896 (2018) argues insertion alone is insufficient.

## Prediction methods

- **IFUM** (bioRxiv 2025.02.10.637420 → Nat Commun 2026, s41467-026-68637-4). Learned model
  (ProtT5 + ESM-IF1 + ESMFold distograms, Evoformer-like), trained on Mega-scale with MMseqs2
  clustering at 0.30 identity. Reports **PCC 0.81 / 0.80 / 0.63** for point mutants / indels /
  double mutants on an 86-domain held-out set, and 0.76 on a 27-domain "common" subset vs ESMtherm
  0.66, ESM2 0.35, FoldX 0.22, Rosetta 0.01. **No statement in the text that this is per-domain**;
  Fig 3a is a single scatter per mutation type, so it is almost certainly pooled.
- **The 0.80 vs 0.62 gap is probably a metric artefact, not a capability gap.** This is the same
  pooled-versus-per-protein trap that killed Gate 9 in this project.
- **pLDDT precedent**: yeast DHFR folding sensor (bioRxiv 2024.07.11.603017, Structure 2024) used
  **delta-pLDDT** (AF2 run on each mutant, minus wild type), one protein, AUC 0.75 insertions /
  0.82 deletions, beating Rosetta ddG, CADD, PROVEAN, ESM1b. Woods et al., Structure 2023
  (PMID 37119820) combined mean pLDDT with Rosetta for deletions.
  **No multi-domain benchmark of wild-type pLDDT at the site, and no comparison with ProDomino,
  was found.** Our 0.845 over 355 domains appears to be unpublished.
- ProteinGym has an indel track (~300k mutants, 74 DMS assays). ESM-1v/ESM-2 masked-marginals
  cannot score indels; the track is autoregressive and family models. Reported leaders
  (*snippet*): ProFam ensemble 0.530, PoET 0.517, DPLM-Evo 0.495, ProGen2-M 0.464, Tranception
  0.463.
- IndeLLM (Patterns 2025/2026): PLL over shared residues, 7,500 clinical indels, MCC 0.77. Clinical
  pathogenicity, not stability.

## Short peptide tags — a real gap, but no data

No library-scale scan of internal GFP11, HiBiT, HA/FLAG/myc, tetracysteine or sortase tags exists
in any protein. Only case studies (tubulin GFP11 in the H1-S2 loop, AF2-guided, PLOS Biol 2024;
EpicTope V5 sites in zebrafish; Abedi 1998, 10 sites in GFP). Pentapeptide/Mu transposon scanning
(Hallet & Sherratt 1997; TRIAD, Nat Commun 2020) exists for a handful of proteins with no public
per-site tables found. **No published correlation between short-insert and domain-insert tolerance
at matched sites.** The best validated site-choice rule is Oesterle et al. 2017 (BMC Biol), which
finds length-variable regions in homolog alignments and filters by secondary structure and
accessibility.

Conclusion: this direction is genuinely open but cannot be tested without new wet-lab data, so it
is a product idea rather than the next experiment.

## What is left unowned

1. The Vetter accommodation magnitudes have never been tested outside T4 lysozyme.
2. Insertion-vs-deletion asymmetry at matched positions, and its dependence on local structure, is
   unmeasured in the Tsuboyama data (Topolska's r = 0.484 is from their own yeast data).
3. The evolutionary and experimental asymmetries point in opposite directions and nobody has
   reconciled them.
4. No one has audited pooled-versus-per-domain metrics across indel predictors, though the field is
   actively comparing numbers computed both ways.
