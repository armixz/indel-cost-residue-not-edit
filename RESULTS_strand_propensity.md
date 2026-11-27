# Does the strand arm have its own residue ordering? No — and that is now a result

Pre-registered in `PREREG_strand_propensity.md`, including the addendum recording the
context precondition before any number was computed. Analysis:
`analysis/51_strand_propensity.py`. Output: `results/strand_propensity.json`.

**Headline: the strong claim does not land. The double dissociation required by the
pre-registration failed on one of its two halves. But the test was worth running, because
it replaces an admission in the manuscript with a measurement, and it adds a specificity
control the paper did not have.**

## Verification, done twice

Both GB1 scales were read from photographic reproductions of the printed journal tables by
the search agent and then **re-read independently by me** from the same scans. All 20 rows
of Minor & Kim Table 1 and all 20 rows of Smith/Withka/Regan Table 2 match. The scans are
kept locally in `refs_scans/`, which is gitignored: they are copyrighted journal pages from
university course handouts and are not redistributed.

The sign trap was real. Minor & Kim and Smith & Regan measured **the same guest position
in the same host protein** and published **opposite sign conventions** — larger = better
for the first, more negative = better for the second. Both were aligned to a penalty scale
(larger = worse) before use, and the pre-registered sanity check that glycine must carry
the highest penalty on every scale passed.

## S5 first, because it licenses everything else

| scale pair | Spearman |
|---|---|
| Pace helix vs Minor & Kim beta | **0.104** |
| Pace helix vs Smith & Regan beta | −0.061 |
| Pace helix vs Kim & Berg beta | −0.026 |
| Minor & Kim vs Smith & Regan | 0.921 |

**The helix and beta scales are essentially uncorrelated.** This matters more than any
other number here: the pre-registration warned that near-collinear scales would leave the
discriminating test underpowered and that any result would then have to be reported as
such. That escape clause does not apply. The test had real power to separate the two
accounts, and what it found is what it found.

The two GB1 beta scales agree with each other at 0.921, as expected from the same host and
position.

## The marks

Correlations of the per-residue contrast against each scale, 18 residues (proline excluded
throughout), paired within domain and averaged over domains:

| contrast | Pace helix | Minor & Kim | Smith & Regan | Kim & Berg |
|---|---|---|---|---|
| strand minus loop | **−0.522** | −0.362 | −0.094 | +0.094 |
| helix minus loop | **−0.639** | −0.226 | −0.057 | −0.004 |
| strand, central only | −0.439 | −0.203 | +0.005 | +0.325 |
| strand, edge only | −0.308 | −0.492 | −0.306 | −0.168 |

**S2 (the test that would have upgraded the claim): FAILED.** On the strand arm, beta
propensity minus helix propensity = **+0.183 [−0.043, 0.411]**, which contains zero. The
point estimate is positive, meaning helix propensity is if anything the *better* predictor
of the strand contrast. There is no evidence that the strand arm has a beta-specific
ordering.

**S3 (the mirror): PASSED.** On the helix arm, helix propensity beats beta propensity by
**−0.374 [−0.586, −0.141]**, excluding zero. The matched scale wins decisively on the
helix arm.

**No double dissociation.** The pre-registration states that one half alone is not enough
and will not be reported as if it were. It is not.

**S4 (scale stability): FAILED, and this is important.** The two beta scales correlate with
each other at 0.921 yet give −0.362 and −0.094 on the same strand contrasts. Whatever weak
beta signal exists is not stable across two scales measured at the same position in the
same protein. This is independent evidence that these scales are not a reliable
instrument, consistent with what the literature already says about them.

**S6 (edge versus central): FAILED TO REPLICATE, in the opposite direction.** The
pre-registered prediction, from Minor & Kim's own context result, was that a propensity
measured at a central strand should work at central strands and fail at edge strands. The
observed difference is **+0.291 [−0.063, 0.669]** — the wrong sign, and containing zero.
The point estimates run the other way: beta propensity correlates at −0.492 at edge
strands and −0.203 at central ones. Reported as a failure to replicate and not
reinterpreted. Note the power asymmetry: 5,137 edge positions against 1,824 central ones.

## What the manuscript gains, and what it does not

**Does not gain:** the stronger claim that each structure class carries its own residue
ordering matched to its own propensity scale. That claim is not supported and will not be
made.

**Gains, first:** an admission becomes a measurement. The manuscript currently says *"No
experimental beta-sheet propensity scale was used, because none could be verified against
its source; the strand contrasts are therefore reported without a propensity analysis."*
Two scales are now verified and tested, and the answer is that they do not predict the
strand contrast better than the helix scale does.

**Gains, second — the real one:** a specificity control. The manuscript reads the helix
propensity scale as a general proxy for contribution to regular secondary structure. The
obvious objection is that any propensity-like scale would do as well, i.e. that the result
is about residue physics in general rather than about propensity. S3 answers it: on the
helix arm, helix propensity beats beta propensity with an interval excluding zero, and the
two scales are uncorrelated (ρ = 0.10), so this is not an artefact of the scales measuring
the same thing. The scale is doing specific work.

**Gains, third:** a sharper statement of the puzzle. Helix propensity predicts both the
helix arm (−0.64) and the strand arm (−0.52). Beta propensity predicts neither well, and
the two families of scale are uncorrelated. So whatever the helix propensity scale is
tracking in this data, it is not beta-sheet-forming ability, and it is not simply
"regular secondary structure" either, or the beta scales would carry some of it. The
honest reading is that the helix scale captures a residue property — plausibly backbone
conformational freedom, with glycine at one extreme — that governs the cost of removing a
residue from any regular element. The manuscript should say that instead of the vaguer
current wording.

## Caveat that constrains all of the above

The context precondition was triggered before any computation and stands: beta-sheet
propensity is strongly context-dependent (Minor & Kim, Nature 1994;371:264, where moving
the guest site within one protein collapses the correlation with statistical preferences
from r = 0.79 to 0.25). A null result against these scales is therefore weak evidence
about the strand arm, not strong evidence. S4's instability confirms the concern
empirically. The correct conclusion is "beta propensity scales do not predict the strand
contrast", not "the strand arm has no residue ordering". The strand arm clearly does have
one — helix propensity predicts it at −0.52 — we simply do not have a sheet-specific
instrument good enough to ask whether a sheet-specific ordering sits underneath.
