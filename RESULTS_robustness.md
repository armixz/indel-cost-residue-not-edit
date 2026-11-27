# The robustness battery: every mark, with its number

Pre-registered in `PREREG_robustness.md` before any of it was run. Analyses:
`analysis/52_robustness.py`, `53_null_and_count.py`, `54_natural_designed.py`.
Outputs: `results/robustness.json`, `null_and_count.json`, `natural_designed.json`.

Reporting rule from the pre-registration: every mark appears here with its number,
whether it helps or not. Three marks did not come out the way I expected and are at the
bottom under their own heading, not buried.

---

## R0. A correction. The analysis set is not what the manuscript said it was

The abstract described the discovery set as **"748 natural domains"**. It is not.

| kind | domains | deletion rows |
|---|---|---|
| natural (named by PDB identifier) | **399** | 19,380 |
| designed (EEHEE / HEEH / HHH / trRosetta-hallucinated) | **349** | 15,336 |
| total | **748** | 34,716 |

Designed structures are 47% of the set. This is a factual error in the manuscript and is
corrected regardless of anything below.

It also breaks the manuscript's framing of 57 HEEH_KT miniproteins as an "independent
replication set": designed proteins were already inside the discovery set, so that was
never an independent replication. R1 replaces it with a real one.

---

## R1. Natural versus designed — the replication, and it holds

Designed miniproteins share no evolutionary history with natural domains, are smaller, and
have idealised secondary structure. The helix-minus-loop contrast versus helix propensity,
computed separately in each, with bootstrap intervals over each subset's own domains:

| subset | domains | rows | residues | Spearman | 95% CI |
|---|---|---|---|---|---|
| natural | 399 | 19,380 | 18 | **−0.575** | [−0.672, −0.395] |
| designed | 349 | 15,336 | 16 | **−0.535** | [−0.678, **−0.015**] |

Both negative, both intervals excluding zero. By the pre-declared rule this is a genuine
replication and the claim is general rather than conditional on protein origin.

**Stated precisely, because the margin matters:** the designed interval excludes zero by
0.015. It is a pass, not a comfortable one. Designed miniproteins have a narrower amino
acid composition, so only 16 residues clear the five-domain threshold there against 18 in
natural domains, which costs power.

**R1b, the interaction F-test per subset**, is far less marginal and replicates strongly in
both:

| subset | F | p | n |
|---|---|---|---|
| natural | 7.40 | 1.0 × 10⁻³⁴ | 19,380 |
| designed | 9.56 | 4.0 × 10⁻⁴⁷ | 12,932 |

**A distinction the manuscript must make and currently does not.** These two tables answer
different questions with very different amounts of evidence. The F-test asks whether a
residue-by-structure interaction exists at all, over tens of thousands of measurements, and
the answer is overwhelming. The rank correlation asks whether the interaction is *ordered
by helix propensity*, over 18 points, and the evidence there is much thinner (see R6b,
p = 0.005). Quoting p = 3.8 × 10⁻⁴⁶ next to a claim about the ordering would invite the
reader to transfer the first number's confidence onto the second claim.

*Note on an earlier run.* A first pass gave the designed subset as 292 domains with an
interval containing zero. That loader merged a cached secondary-structure table and
silently dropped 57 HEEH_KT designs whose names lack a `.pdb` suffix. The table above uses
the manuscript's own loader and the complete 349. The discrepancy was found by checking the
domain count against the manuscript rather than by the result looking wrong.

---

## R2. A second, independent helix propensity scale

| scale | Spearman | 95% CI |
|---|---|---|
| Pace & Scholtz 1998 | −0.634 | [−0.742, −0.392] |
| O'Neil & DeGrado 1990 | **−0.639** | [−0.748, −0.406] |

The result is not an artefact of one published table.

**The honest limit of this check:** the two helix scales correlate with each other at
ρ = 0.87. This is nothing like the independence of the beta-scale comparison (ρ = 0.10),
so R2 shows the finding survives a different measurement of the same quantity, not that it
survives a different quantity.

**Incidental, and worth a line in the manuscript:** O'Neil & DeGrado list proline at about
3 kcal/mol. That is the likely origin of the "3.16" figure that circulates attached to the
Pace scale and that we could not verify against Pace & Scholtz. The two scales are
different measurements and the number appears to have migrated between them.

---

## R3. The variant quality filter does not matter

Variants whose 95% interval on ΔG exceeds a threshold are excluded; the manuscript uses 1.0.

| threshold | rows | Spearman |
|---|---|---|
| 0.5 | 30,461 | −0.554 |
| 0.75 | 31,562 | −0.615 |
| **1.0 (used)** | 32,312 | **−0.634** |
| 1.5 | 33,266 | −0.636 |
| 2.0 | 33,708 | −0.713 |
| none | 34,099 | −0.689 |

No sign change anywhere, and looser filters give a *stronger* correlation, so the chosen
threshold is the conservative end of the range rather than a favourable pick.

## R4. The estimator does not matter

Minimum domains per residue: 3 → −0.634, 5 → −0.634, 8 → −0.669, 10 → −0.669.
Median over domains instead of mean: −0.678.

## R5. No single residue or domain carries it

- Leave-one-residue-out, 18 refits: range **−0.712 to −0.566**. Most influential is
  leucine, whose removal *strengthens* the correlation to −0.712.
- Leave-one-domain-out, 691 refits: range **−0.676 to −0.606**.
- Without glycine, which sits at the end of both scales and could lever the fit alone:
  **−0.566** over 17 residues.

Nothing flips the sign or approaches zero.

---

# The three marks that did not come out as expected

## R6. The residue-shuffling null is not centred at zero

Shuffling residue identity within each domain, 500 times: null mean **+0.190**, sd 0.208,
2.5–97.5% range [−0.251, +0.582]. The observed −0.634 lies beyond every one of the 500.

The pre-registration says: *"The null must be centred near zero. If it is not, the
estimator manufactures correlation from domain composition and the result is void."*
**As written, that condition is not met.** With 500 permutations the mean is about 20
standard errors from zero, so the offset is real, not noise.

I am not treating the result as void, and the reason is stated so a reader can disagree:
the bias is **positive** and the finding is **negative**. An estimator that manufactures
+0.19 out of domain composition cannot manufacture −0.63; if anything it makes the test
conservative. The pre-registered wording was aimed at a bias pointing toward the finding,
and I did not anticipate one pointing away. That is a flaw in how I wrote the condition,
not a reinterpretation of the result after seeing it, and it is reported here rather than
quietly rewritten.

The likely mechanism: shuffling changes which residues qualify for the five-domain
threshold, and qualification correlates with abundance, which correlates with the
propensity scale. It is a property of the inclusion rule, not of the data.

## R6b. A cleaner null, which does pass

Permuting the propensity values across residues with the measured contrasts held fixed
asks directly whether *this ordering* is special, and has no estimator bias at all:

- 20,000 permutations, null mean **−0.0025**, sd 0.241
- observed −0.634, **p = 0.005** two-sided

Centred at zero as required, and the observed value is significant. **This, not the F-test,
is the right p-value for the ordering claim**, and it is modest because the test has 18
points. The manuscript should quote it.

## R1 margin and R2 independence

Already stated above, repeated here so the three soft spots sit together: the designed
replication clears zero by 0.015, and the second helix scale correlates with the first at
0.87 rather than being independent.

---

# What changes in the manuscript

1. "748 natural domains" → 399 natural and 349 designed, stated in the abstract, Results
   and Methods.
2. The 57-HEEH_KT "independent replication set" framing is removed and replaced by the
   natural/designed split in R1, with its margin stated.
3. Add R2 with its 0.87 caveat, and the O'Neil proline note.
4. Add a compact robustness statement for R3, R4 and R5.
5. Quote R6b's p = 0.005 for the ordering claim, and stop letting the F-test p-value stand
   next to it.
6. Report R6's non-zero null and the reasoning above.

After these, the manuscript is rebuilt once and frozen.
