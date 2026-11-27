# Pre-registration: does the strand arm have its own residue ordering?

Written 2026-11-26, BEFORE any beta-sheet propensity value was obtained and before any
test was run. The literature search for the scales was running while this was written;
nothing from it has been read.

## The gap this addresses

The manuscript reports that the extra cost of deleting a residue from a helix rather than
a loop tracks that residue's experimental helix propensity. It then hedges, twice:

1. Helix propensity predicts the **strand**-minus-loop contrast nearly as well, and the
   helix-minus-strand difference does not correlate reliably with it. We therefore read
   the scale as a general proxy for "contributes to regular secondary structure" rather
   than as a helix-specific mechanism.
2. No experimental beta-sheet propensity scale was used, because none could be verified
   against its source. The strand contrasts are reported with no propensity analysis.

Both hedges exist because one axis is missing. If the strand contrast tracks **beta**
propensity better than it tracks **helix** propensity, the claim strengthens from "residue
identity matters" to "each structure class has its own residue ordering, matched to its
own propensity scale". If it does not, the current general-proxy reading is confirmed and
the manuscript is unchanged.

## Precondition, which may stop this before it starts

A scale is usable here only if its values are verified against the primary source, with
the sign convention quoted from the paper itself. Values from AAindex, review tables or
other secondary compilations are not usable, for the same reason the manuscript already
refuses them. If no scale clears that bar, this test does not run and the manuscript keeps
its current wording. That outcome is declared acceptable in advance.

A second precondition: beta-sheet propensity is known to be strongly context-dependent
(Minor & Kim's own follow-up is titled "Context is a major determinant of beta-sheet
propensity"). If the published context effect is comparable in size to the spread of the
scale itself, then a single context-free beta scale is not a meaningful predictor, and
a null result here would be uninformative rather than evidence against the hypothesis.
Whether that is so will be determined from the literature BEFORE the marks are read, and
recorded here, so that a null cannot be reinterpreted afterwards.

---

# Addendum written after the literature search, before any mark was computed

## The context precondition is TRIGGERED. Recorded here before any number was read.

The literature search is complete and no analysis has been run. Applying the precondition
above, honestly:

**Minor & Kim 1994b (Nature 371:264, PMID 8078589) is decisive against a context-free
beta scale.** In the same host protein, same lab, same method, moving the guest site from
the central strand to the edge strand destroys the relationship: correlation with
statistical beta preferences falls from r = 0.79 (central, position 53) to r = 0.25
(edge, position 44), as quantified by Avbelj & Baldwin (PNAS 2002;99:1309). Their own
conclusion, verbatim from the abstract: *"beta-sheet formation is determined in large part
by tertiary context, even at solvent-accessible sites, and not by intrinsic secondary
structure preferences."*

Corroborating: Otzen & Fersht (Biochemistry 1995;34:5718) conclude *"thermodynamic scales
for beta-sheet propensity are currently of insufficient precision for general design
purposes."* The Ala-to-best-residue gap differs by a factor of 5-8 between host systems
(0.21 kcal/mol in the zinc finger, 1.63 in GB1), i.e. the between-system disagreement is
as large as the within-scale span.

**Consequence, fixed now:** a null result on S1 or S2 is UNINFORMATIVE. It cannot be
reported as evidence that the strand arm lacks a residue ordering, because the predictor
is known to be unreliable. Only a POSITIVE result is interpretable, and it would be
interpretable precisely because it was obtained despite a noisy predictor.

## Scales admitted, and the sign trap

Verified against photographic reproductions of the printed tables:

| scale | host, position | reference | convention as printed |
|---|---|---|---|
| Minor & Kim 1994a, Nature 367:660 | GB1, residue 53, **central** strand | Ala = 0 | **larger = better** beta former |
| Smith, Withka & Regan 1994, Biochemistry 33:5510 | GB1, residue 53, **central** strand | Ala = 0 | **more negative = better** |
| Kim & Berg 1993, Nature 362:267 | zinc finger, solvent-exposed antiparallel sheet | Gly = 0 | **more negative = better** |

Kim & Berg is verified at one remove (reproduced in Smith et al.'s own Table 2 and matched
independently), and its dynamic range is 4-5x compressed relative to the GB1 scales; it is
reported as a third check, never as the primary.

**Minor & Kim and Smith & Regan use OPPOSITE sign conventions while measuring the same
guest position in the same host.** To align everything with the manuscript's helix scale,
which is a *penalty* where larger = worse former:

- `beta_penalty(Minor-Kim) = -(printed value)`
- `beta_penalty(Smith-Regan) = +(printed value)`
- `beta_penalty(Kim-Berg) = +(printed value)`

A sanity check is fixed in advance: after alignment, glycine must carry a high penalty on
all three scales. If it does not, the alignment is wrong and the analysis stops.

Proline is excluded throughout: Smith & Regan report it as not determined, and Minor & Kim
give only a bound.

**Expected sign.** On the helix arm the contrast correlates NEGATIVELY with the Pace
penalty: a poor helix former contributes little to the helix, so deleting it from one
costs little extra. The strand arm should behave the same way, so S1 predicts a NEGATIVE
correlation between the strand-minus-loop contrast and the aligned beta penalty.

## S6, added before any computation: the edge-versus-central test

The context finding is not only a limitation; it is a prediction we can test, and it is
the reason this analysis is now worth running even though a null is uninformative.

Both admitted GB1 scales were measured at a **central** strand. Our strand positions are a
mixture of central and edge strands. Minor & Kim's result says a propensity measured at a
central strand should apply at central strands and fail at edge strands. So:

- **S6.** Split our strand positions by bridge topology from `analysis/17_dssp.py`: a
  residue paired with partners on both sides is central, one side only is edge. Recompute
  S1 separately in each class.
  - **Predicted:** the beta propensity correlation is stronger at central strands than at
    edge strands.
  - If that holds, it reproduces a 1994 result obtained on two guest positions of one
    protein across hundreds of domains and tens of thousands of deletions, in a completely
    different kind of experiment. That is a finding in its own right, independent of
    whether S2's double dissociation succeeds.
  - If it does not hold, it is reported as a failure to replicate, and it is not
    reinterpreted as support for a context-free scale.

S6 is declared before the data are touched, and its prediction is fixed in the direction
stated, so it cannot be flipped afterwards.

## Estimator

Identical to the manuscript's, with helix replaced by strand. For each domain and each
amino acid, the mean deletion cost at strand positions minus the mean at loop positions;
averaged over domains; bootstrap intervals resampling domains. Proline excluded, as
throughout. Secondary structure from analysis/17_dssp.py.

## Marks

- **S1 (the test).** Spearman correlation of the strand-minus-loop contrast against
  verified beta-sheet propensity.
  - The sign expected is the same as the helix arm: a residue that is a poorer beta-sheet
    former should cost less to delete from a strand, so the contrast should fall as the
    propensity penalty rises. Whether that means a negative or positive rho depends on
    each scale's sign convention, which is why the convention must be quoted from the
    source before the number is interpreted.

- **S2 (the comparison that decides it).** The difference between two correlations on the
  SAME strand contrasts: against beta propensity, and against helix propensity. Bootstrap
  the difference over domains.
  - **Strengthens the paper** only if the beta correlation is larger in magnitude AND the
    bootstrap interval on the difference excludes zero.
  - **Confirms the current hedge** if the interval contains zero, i.e. the two scales
    predict the strand contrast equally well. This is the outcome the manuscript already
    assumes, and it requires no change.
  - **Weakens the paper** if helix propensity beats beta propensity on the strand arm with
    an interval excluding zero. That would mean the scale is tracking something other than
    class-specific propensity, and the Results paragraph would need rewriting.

- **S3 (the mirror).** The same comparison on the HELIX arm: does helix propensity beat
  beta propensity there? The strong claim requires both halves, a double dissociation.
  One half alone is not enough and will not be reported as if it were.

- **S4 (scale agreement control).** If more than one verified beta scale exists, run S1
  with each. A result that holds under one scale and not another is reported as unstable,
  not as a finding under the scale that worked.

- **S5 (collinearity check).** Spearman between the beta and helix propensity scales
  themselves. If the two scales are strongly correlated, S2 has little power to separate
  them, and that must be stated with the result rather than discovered later. Computed and
  reported regardless of outcome.

## Pre-declared interpretation

- S2 and S3 both favour the matched scale, both intervals excluding zero: double
  dissociation. The manuscript's claim upgrades and the Results paragraph is rewritten.
- S2 contains zero: the general-proxy reading stands. The manuscript gains one verified
  sentence saying the strand arm was tested against a beta scale and the two scales do not
  separate, which closes the current "no scale could be verified" admission with a result
  instead of an absence.
- S5 shows the scales are near-collinear: whatever S2 returns is reported as
  underpowered, and no dissociation claim is made in either direction.

## What will not be done

No scale will be selected after seeing which one gives the better correlation. The scales
to be used are fixed by the verification step alone, before any mark is computed, and
every verified scale that clears the precondition is reported in S4 whether it helps or
not.
