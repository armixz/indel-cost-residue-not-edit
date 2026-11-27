# Pre-registration: how a helix absorbs an indel, at mega-scale

Written 26 Nov 2026, **before computing anything**. Literature reconnaissance (five searches)
was completed first and is summarised in LITERATURE_indel.md. Nothing below was chosen after
seeing a result in this dataset.

## Why this question

Vetter, Baase, Heinz, Xiong, Snow & Matthews (Protein Sci 1996, 5:2399, PMID 8976549) determined
how an alpha-helix accommodates an inserted residue, from 32 insertions and 5 deletions in nine
helices of T4 lysozyme, with structures for eight variants. Their abstract gives three class
averages of destabilization:

| accommodation mode | destabilization |
|---|---|
| translocation (register shift, spilled into the adjacent loop) | 2.5 kcal/mol |
| looping-out near a helix end | 2.6 kcal/mol |
| bulging within the body of the helix | 4.9 kcal/mol |

Heinz, Baase & Matthews (J Mol Biol 1994, 236:869) state the mechanism: "The inserted amino acids
were incorporated into the helix by replacing preceding wild-type amino acids and causing a shift
in register towards the N terminus. As a consequence, wild-type amino acids were translocated from
the helix into the preceding loop."

**These numbers have never been tested outside T4 lysozyme.** Tsuboyama et al. 2023 (Zenodo
7992926) Lib3 contains 95,184 single-residue insertion and 46,660 single-residue deletion
measurements across ~400 domains, in the same unit (kcal/mol). Topolska, Beltran & Lehner
(Nat Commun 2025, 16:2617) analysed the same insertions and report, qualitatively, that "both ends
of alpha-helices are indeed slightly more tolerant of insertions, deletions and substitutions",
and that bulge accommodation "is a rare occurrence". They give **no kcal/mol value** for the
end-versus-middle difference, and they do not test the Vetter model's magnitudes.

A live contradiction also exists. Yang et al. 2024 (PMC11102076, 11,444 alignments) find from
evolution that **insertions** are the more depleted class inside secondary structure (O/E in
strands 0.26-0.39 for insertions versus 0.41-0.52 for deletions), while Topolska find from
experiment that **deletions** are the more damaging class overall (68.8% versus 62.1% strongly
reduced). Nobody has reconciled these. This dataset carries insertions and deletions at matched
positions, which is what a reconciliation requires.

## Sign convention

`deltaG` in Lib3 is folding stability, positive = folded (median wild-type 2.52 kcal/mol).
ddG = deltaG(variant) - deltaG(wild type). **Negative ddG = destabilizing.** Vetter's
"destabilization of 4.9 kcal/mol" is ddG = -4.9 here. Below, "cost" means -ddG, so larger cost =
worse.

## Prerequisite (not a test)

Secondary structure in this project has so far been assigned from phi/psi regions only. That is
too coarse here, because the question is about helix boundaries. Before any mark is evaluated I
will implement Kabsch-Sander hydrogen-bond secondary structure (DSSP rules for H and E) on the
Tsuboyama AlphaFold models, define helix elements as maximal runs of H of length >= 5, and record
for every position: element type, index within the element, distance to the nearest element
terminus, and distance to the nearest loop residue. The assignment will be sanity-checked against
the existing phi/psi labels (expected agreement >80%; if below, I stop and investigate rather than
proceed).

---

## H1. Does the Vetter magnitude hold outside T4 lysozyme?

Vetter: a bulge in the helix body costs 4.9, looping-out near a helix end costs 2.6. If accommodation
mode is set by position within the helix, insertion cost should be **about 2.3 kcal/mol higher in
the middle of a helix than near its ends**.

Define: helix END = index 1 or 2 from either terminus of the element; helix MIDDLE = index >= 4
from both termini. Statistic: mean cost(MIDDLE) - mean cost(END) for single-Ala insertions in
helices, with a 95% CI from bootstrapping **over domains** (never over positions).

- **Mark H1a (direction and a non-trivial size).** Difference >= 1.0 kcal/mol and the CI lower
  bound > 0.5. This is deliberately set below Vetter's 2.3, so passing it does not require the
  crystallographic numbers to transfer exactly.
- **Mark H1b (the Vetter magnitude itself).** The 95% CI contains 2.3 kcal/mol.
- **Mark H1c (not a burial artefact).** H1a's difference holds, same sign, in every tertile of
  relative solvent accessibility, and in every tertile of pLDDT. Helix middles are more buried than
  helix ends, so without this control the result is uninterpretable. **If H1c fails, H1 fails**,
  whatever H1a did.

Expected tension, stated now: Topolska's word "slightly" suggests the real difference may be well
below 1.0, in which case H1a fails and the honest finding is that the T4 magnitudes do **not**
generalise. That is a result, and it will be reported as plainly as a pass.

## H2. Insertion and deletion are absorbed by different resources

Mechanism from Heinz/Vetter: an insertion must put extra residues somewhere, and they are
translocated into the adjacent loop; a deletion is absorbed by a loop at the helix end being
shortened. Both need an accessible loop, but only the insertion must **deposit** material there.
So loop proximity should help insertions more than deletions.

Statistic: per position inside helices, asymmetry A = cost(insertion) - cost(deletion), using the
Ala insertion. Positive A = insertion worse.

- **Mark H2a.** A decreases with proximity to the nearest loop: Spearman rho between A and distance
  to nearest loop >= 0.15 in magnitude, with a by-domain bootstrap CI excluding 0.
- **Mark H2b (sign flip).** Mean A < 0 at helix positions within 2 residues of a loop, and mean A > 0
  at positions >= 5 residues from any loop, with non-overlapping by-domain CIs.
- **Mark H2c (residue-identity control).** A deletion removes the wild-type residue, while the
  insertion is always Ala, so A confounds "where" with "what". H2a must survive in the subset where
  the wild-type residue **is** alanine (a true like-for-like comparison), with the same sign. If
  the alanine-only subset is too small for a CI (fewer than 200 positions or 30 domains), H2 is
  reported as **untested**, not as passed.

## H3. Reconciling evolution with the bench (exploratory, no pass/fail)

Compute A by secondary structure class (helix, strand, loop). Yang et al. find insertions more
depleted than deletions inside structure, most strongly in strands; Topolska find deletions more
damaging overall. These are compatible only if the asymmetry reverses by context. I will report A
per class with by-domain CIs and state whether it reverses. **This is labelled exploratory in
advance** and cannot be reported as a confirmed finding from this dataset alone.

## B. Metric audit (separate, near-certain to yield a number)

INDELi (Topolska) reports per-domain mode r of 0.620 (deletions) and 0.608 (insertions) over 171
domains. IFUM (Nat Commun 2026, s41467-026-68637-4) reports Pearson 0.80 for indels on a held-out
Mega-scale set, from what appears to be a single pooled scatter over 86 domains. A pooled
correlation absorbs between-domain variance in mean ddG and is not comparable to a per-domain one.
This project was itself misled by exactly this difference in Gate 9.

- **Mark B1.** For the same predictor and the same indel variants, pooled Pearson r exceeds the
  per-domain median r by >= 0.10. Predictors used: pLDDT, relative solvent accessibility, and
  distance to terminus (all already computed here; no new model is trained).
- B2 (descriptive). Report the full per-domain distribution, including the lower decile, for each.

B does not depend on H1-H3 and will be reported even if every other mark fails.

## Negative controls, fixed now

1. **Permutation.** Shuffle ddG within each domain and re-run H1 and H2. Every effect must vanish
   (CI containing 0). If any survives, the pipeline is wrong and nothing is reported.
2. **Gly arm.** H1 and H2 are pre-specified on Ala insertions. The Gly insertions are an independent
   replication arm, run **after** the Ala arm is complete and reported either way. Gly has different
   helix propensity, so a difference between arms is informative, not a failure.
3. **Domain-level reporting.** Every CI is bootstrapped over domains. No statistic is reported
   pooled over positions without its per-domain counterpart beside it.

## What would make me stop

- Agreement between the new DSSP assignment and the existing phi/psi labels below 80%.
- Any permutation control that does not vanish.
- Fewer than 30 domains contributing to a given mark.

## Order of execution

1. DSSP-quality annotation + sanity check (prerequisite).
2. B (metric audit) — independent of the mechanism.
3. H1 with its H1c control.
4. H2 with its H2c control.
5. Permutation controls.
6. H3 exploratory, and the Gly replication arm.

Results go to RESULTS_indel.md, marks reported as passed or failed exactly as written above.
