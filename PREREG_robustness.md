# Pre-registration: the full robustness battery, run once before submission

Written 2026-11-26, before any of the tests below were run. Nothing here is reported
selectively: every mark defined in this file appears in the results write-up with its
number, whether it helps or not.

## Why this exists

The manuscript has been edited several times in response to single findings. That is a bad
way to finish a paper: each edit is defensible on its own, but the sequence looks like
tuning, and the author cannot tell which choices the result actually depends on. This
document fixes the complete list of remaining checks in advance, so they are run as one
batch and reported as one batch.

After this battery, the manuscript is edited once and then frozen.

## R0. A correction that is not a test

The analysis set of 748 domains is described in the abstract as "748 natural domains".
That is wrong. It is **399 natural domains** (named by PDB identifier), **252 designed
miniproteins** (EEHEE, HEEH, HHH and similar) and **97 trRosetta-hallucinated designs**
(`*_TrROS_Hall`), i.e. 349 designed structures in total, 47% of the set.

This must be corrected in the abstract, the Results and the Methods regardless of what
the tests below show. It also means the manuscript's current framing of 57 HEEH_KT designs
as an "independent replication set" is misleading, because designed proteins are already
inside the discovery set. The replication framing has to be rebuilt around the split in R1.

## R1. Natural versus designed — the replication test that matters

Designed miniproteins share no evolutionary history with natural domains, are smaller, and
have idealised secondary structure. If the residue-by-structure interaction is a real
biophysical effect it should appear in both. If it appears only in one, the manuscript's
central claim is conditional and must say so.

- **R1a.** Helix-minus-loop contrast versus helix propensity, computed separately in the
  natural subset and the designed subset, each with bootstrap intervals over its own
  domains.
- **R1b.** The interaction F-test, separately in each subset.
- **R1c.** The deletion-minus-insertion offset, separately in each subset.

**Pre-declared interpretation.** Both subsets negative with intervals excluding zero:
the claim is general and the manuscript gains a genuine replication. One subset only:
the claim is conditional on protein origin and the abstract must be rewritten to say which.
Opposite signs: the effect is an artefact of pooling two populations and the paper does
not go out.

## R2. A second, independent helix propensity scale

The manuscript rests on one scale (Pace & Scholtz 1998). A reviewer will ask whether the
result survives a different one. O'Neil & DeGrado (Science 1990;250:646) is an independent
experimental helix propensity scale measured in a different host by a different method; its
printed table was read from the same scan that carried the Smith & Regan table, and its
convention is already "larger = worse former" (Ala $-0.77$ best, Gly $0.00$, Pro $\approx 3$).

- **R2a.** Helix-minus-loop contrast versus the O'Neil--DeGrado scale.
- **R2b.** Spearman between the two helix scales, so the reader can judge how independent
  the check really is.

**Pre-declared interpretation.** If R2a is negative with an interval excluding zero, the
result is not an artefact of one scale. If it is null, the manuscript must say the result
depends on the choice of helix scale, which would be a serious weakness.

Incidental: O'Neil & DeGrado list proline at about 3 kcal/mol. This is the likely origin of
the "3.16" value that circulates for the Pace scale and that we could not verify. The
manuscript should say so, because it explains a discrepancy a reader may otherwise trip over.

## R3. The variant quality filter

Variants whose reported 95% interval on $\Delta G$ exceeds 1.0 kcal/mol are excluded. That
threshold was chosen once and never examined.

- **R3.** Repeat the headline correlation at thresholds 0.5, 0.75, 1.0, 1.5, 2.0 and with
  no filter at all, reporting the number of variants retained at each.

**Pre-declared interpretation.** The conclusion must not depend on the threshold. If the
correlation changes sign or loses its interval anywhere in that range, the threshold is
doing work and the manuscript must report the whole sweep rather than one value.

## R4. The estimator

Two arbitrary choices: a residue enters the correlation only if it has a contrast in at
least 5 domains, and the per-residue contrast is the mean over domains.

- **R4a.** Minimum-domain threshold at 3, 5, 8 and 10.
- **R4b.** Median over domains instead of mean.

## R5. Influence

A rank correlation over 18 points can be moved by one residue or a few domains.

- **R5a.** Leave-one-residue-out: recompute the correlation 18 times, report the range and
  which residue moves it most.
- **R5b.** Leave-one-domain-out over all domains, report the range.
- **R5c.** The correlation with glycine excluded, since glycine sits at the end of both the
  propensity scale and the contrast range and could lever the fit alone.

**Pre-declared interpretation.** If removing any single residue or domain flips the sign or
moves the correlation outside the reported bootstrap interval, the result is fragile and
the manuscript says so.

## R6. A permutation null for this specific correlation

The manuscript has a permutation control for the insertion contrast but not for the
propensity correlation itself.

- **R6.** Shuffle residue identity within each domain, recompute the per-residue contrasts
  and the correlation against helix propensity, 500 times. Report the mean and standard
  deviation of the null and where the observed value sits.

**Pre-declared interpretation.** The null must be centred near zero. If it is not, the
estimator manufactures correlation from domain composition and the result is void.

## Reporting rule

Every mark above appears in `RESULTS_robustness.md` with its number. Marks that fail are
reported in the manuscript, not only in the repository. After this battery the manuscript
is edited once, rebuilt once, and not changed again before external review.
