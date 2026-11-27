# Pre-registration 2: why are deletions worse than insertions in beta strands?

Written 26 Nov 2026, **after** the H3 result and **before** computing anything below.

## Status of the finding so far

PREREG_indel.md declared H3 (asymmetry by structure class) **exploratory**. It gave, on Tsuboyama
Lib3, the per-position asymmetry A = cost(Ala insertion) - cost(deletion), averaged within domain
then over domains:

| class | A (kcal/mol) | 95% CI by domain |
|---|---|---|
| beta strand | **-0.587** | [-0.709, -0.467] |
| alpha helix | -0.266 | [-0.331, -0.200] |
| loop | -0.243 | [-0.296, -0.195] |

Negative = the insertion is cheaper. Controls already run (22_asym_control.py, 19,170 positions,
387 domains) rule out the obvious confounder:

- restricted to positions whose wild-type residue **is alanine**, so the deleted and the inserted
  residue are chemically identical, the strand value becomes **-0.842 [-1.133, -0.561]** - larger,
  not smaller
- present in all three residue-volume tertiles (-0.868 / -0.451 / -0.683), with no volume trend
- strand more negative than non-strand within 13 of 16 wild-type residue types, median -0.35
- OLS with domain-clustered errors: is_E coefficient -0.198 after volume, rSASA and pLDDT
- **burial dependence differs by class**: in helices and loops the asymmetry nearly vanishes at
  exposed positions (-0.086 and -0.153) while in strands it persists (-0.571)

That last line is the reason for this document. In helices and loops the asymmetry looks like
ordinary packing - removing a buried residue costs more than adding a small one. In strands
something survives exposure, so something specific to beta sheets is at work.

Cross-check (23_ss_crosscheck.py): against the Lehner lab's own structure table for these domains,
three-state agreement is 0.801 with **zero** H/E confusion, rSASA correlates at 0.955 and element
length at 0.941. With **their** labels the ordering holds but is smaller (strand -0.568 against
loop -0.416); with the stricter labels used here it is larger (-0.763 against -0.388). The
magnitude is therefore sensitive to how strictly a strand is defined, and both numbers must be
reported together.

## The hypothesis

A **beta bulge** is an extra residue on one strand of a beta ladder, accommodated while the
hydrogen-bond ladder on both sides is preserved. Bulges are common: Chan, Mazzarella & Thornton
1993 found 362 in 170 proteins; a 2021 PDB-wide survey counts 5,133 classic and 4,854 G1 bulges.
There is no corresponding structure that absorbs a **missing** residue: a deletion must close the
gap by pulling the chain, which breaks the register over the rest of the ladder.

**H-bulge: an inserted residue in a beta strand has an accommodation route that a deleted residue
does not, so the insertion/deletion asymmetry is larger in strands than in any other class.**

Note the tension this creates with the literature, which is why it is worth testing: Topolska et
al. state that bulge accommodation "is a rare occurrence" in these same data, and Yang et al. 2024
find from evolution that insertions are the **more** depleted class inside structure, most
strongly in strands (O/E 0.26-0.39 for insertions against 0.41-0.52 for deletions). H-bulge
predicts the opposite of the evolutionary pattern, so if it holds, the reconciliation has to be
that evolutionary depletion and thermodynamic cost are measuring different things.

## Marks

### Replication on an independent assay (run first; if this fails, the rest is reported as a
### single-dataset observation and nothing more)

Topolska et al.'s own yeast abundancePCA data (Zenodo 10830272, `scaled_variants_aPCA.rds`):
1,450 deletions and 1,420 insertions across 9 domains, with their own structure labels, rSASA, and
a per-variant error. The assay is orthogonal to Tsuboyama's in-vitro proteolysis.

- **R1.** A is more negative in strands than in non-strand positions, same sign in at least **7 of
  9** domains, with the mean difference's by-domain bootstrap CI excluding 0.
- **R2.** The strand excess survives restricting to **exposed** positions (rSASA above the median),
  where the helix and loop asymmetry was shown to vanish in Tsuboyama. This is the diagnostic mark:
  it separates "beta sheets are special" from "burial".
- **R3 (identity control).** Same sign when restricted to positions where the inserted residue
  equals the wild-type residue. Reported as untested if fewer than 100 such positions.

### Mechanism, on Tsuboyama Lib3

These need bridge-partner information, which 17_dssp.py computes but does not currently save.
24_bridges.py will record, per strand residue: number of bridge partners (1 = edge strand,
2 = interior), bridge type (parallel / antiparallel), and partner identity.

- **M1 (antiparallel).** Classic beta bulges are overwhelmingly antiparallel (Richardson 1978;
  Chan 1993). If bulge accommodation drives the asymmetry, A is more negative in **antiparallel**
  bridges than in parallel ones. Mark: difference <= -0.15 kcal/mol with a by-domain CI excluding 0.
- **M2 (edge strands).** An edge strand has a free side and can bulge outward with less cost than
  an interior strand hemmed in on both sides. Mark: A more negative at **edge** than at interior
  positions, difference <= -0.15 with CI excluding 0.
- **M3 (position within the strand).** A bulge is a local accommodation; closing a deletion's gap
  propagates along the ladder. So the asymmetry should **grow toward the middle** of a strand.
  Mark: A at strand positions with idx_end >= 2 more negative than at idx_end = 0, CI excluding 0.

M1 and M2 are the load-bearing ones. **If both fail, H-bulge is rejected** and the strand excess is
reported as real but mechanistically unexplained.

### Kill conditions, fixed now

1. **Strand length.** Strands are short (median 5 here). If the strand excess disappears once
   element length is matched between classes, it is a length artefact, not a bulge effect. Run as
   a mandatory control, not an optional one.
2. **Permutation.** Shuffling cost within domain must abolish M1-M3.
3. **Measurement error.** Deletions and insertions must have comparable `dG_ci`; if deletions are
   systematically noisier, the asymmetry could be regression to the mean. Check and report the two
   error distributions. If they differ materially, re-run the headline on a tight-CI subset.
4. **The 2-residue strand floor.** 17_dssp.py demotes strand runs shorter than 2 to loop. Re-run
   the headline with a floor of 3 to confirm the result is not created by the shortest strands.

### What I expect

Honest statement before running: the evolutionary evidence points the other way, and Topolska
looked at bulges in these data and called them rare. I give H-bulge less than even odds. The most
likely outcome is that the strand excess is real (the controls are already strong) but that M1 and
M2 fail, leaving it unexplained. That outcome will be reported as plainly as a success.
