# The secondary-structure dependence of indel cost belongs to the residue, not to the indel

**142,000 single-residue insertions and deletions, compared within proteins**

Draft, 26 November 2026. Numbers are final; framing and authorship are not.
This draft is the third version of the night. Two earlier framings were withdrawn after their own
controls failed; both are kept in §6 because the failures are the most useful part.

---

## Abstract

Single-residue insertions and deletions are the second most common coding variant and the least
well characterised. Using the mega-scale folding-stability measurements of Tsuboyama et al. (2023)
we compare the two edits **at the same position in the same protein**, across 18,699 positions in
385 natural domains and, as an independent replication, 2,402 positions in 57 designed HEEH
miniproteins that appear not to have been analysed for indels before.

**(i) A deletion costs more than an insertion at the same site, by 0.303 kcal/mol [0.257, 0.351]
on average.** The direction is already known qualitatively; the paired quantity is not.

**(ii) The penalty shrinks about twofold from tolerant to intolerant positions, from 0.52 to
0.26 kcal/mol, but it does not vanish.** Measured by binning on an *independent* second insertion
at the same position; binning on the same measurement, as is natural to do, produces a spurious
decline all the way through zero, and an errors-in-variables correction moves the regression slope
of deletion on insertion cost from 0.71 to 0.87.

**(iii) The secondary-structure dependence of indel cost belongs to the residue involved, not to
the indel.** On the deletion side, where all 20 residues occur, residue identity interacts with
structure class at F = 7.4, p = 1.7e-36 over 19,897 deletions in 387 domains with domain-clustered
errors and after burial and volume, and the interaction vanishes when residue identity is permuted
within domain (p = 0.68). The interaction has the predicted shape: the extra cost of deleting a
residue from a helix rather than a loop scales with how good a helix former it is, Spearman
**-0.603** against Pace & Scholtz helix propensity across 19 residues, bootstrap CI over domains
[-0.691, -0.446]. On the insertion side only alanine and glycine exist in these libraries, and
switching between them **reverses the sign** of the helix contrast (-0.082 to +0.260) and removes
most of the strand contrast (-0.324 to -0.066); the same switch abolishes the helix contrast in
the designed miniproteins. A statement of the form "strands tolerate indels worse than helices"
therefore has no determinate value until the residue is specified.

We also show that the crystallographic estimates of helix indel accommodation from T4 lysozyme do
not transfer (0.481 [0.387, 0.575] where 2.3 was predicted), and that pooled and per-domain
correlations currently used to compare indel predictors cannot be ordered against each other.

---

## 1. What was already known

Topolska, Beltran & Lehner (Nat Commun 2025;16:2617) report 1-aa insertions as slightly better
tolerated than deletions (KS p = 0.034) and give separate secondary-structure effect sizes for each
edit (strand against helix: Cohen's d = -0.172 for deletions, -0.117 for alanine insertions). They
report a correlation of r = 0.484 between the two edits in their yeast data, with no slope or
functional form. Nielsen et al. (bioRxiv 2024.07.11.603017), on human DHFR, state that "deletions
were more detrimental in beta-strands than in alpha-helices, which was not the case for
insertions" and suggest it is "easier to accommodate an insertion, for example by extending a
loop or forming a bulge, than to handle the increased strain caused by a deletion".

**We do not claim (i).** What is new is the paired quantification, the tolerance dependence in
(ii), the residue dependence in (iii), the replication set, and three negative results.

## 2. Methods

**Data.** Tsuboyama et al. 2023 (Zenodo 7992926) library 3: 95,184 single-residue insertions
(alanine and glycine) and 46,660 deletions. cost = dG(wild type) - dG(variant); variants with a dG
confidence interval wider than 1 kcal/mol are dropped. 18,699 positions carry all three
measurements.

**Structure.** Kabsch-Sander hydrogen-bond assignment on the deposited AlphaFold models, validated
against geometry it does not use (helix centres at phi -63.0, psi -41.8, CA(i)-CA(i+4) 6.2 A;
strand at -113.7, +134.1, 12.2 A) and against the Lehner lab's independent table for the same
domains (three-state agreement 0.801, zero helix/strand confusion, rSASA r = 0.955).

**Replication set.** 69 designed HEEH_KT miniproteins in the same library carry single indels; 57
have deposited models, ~131 variants each. Numbering verified at 1.00 over 2,412 deletion
variants. We find no published analysis using them.

**Statistics.** Every class comparison is computed inside each protein and then averaged over
proteins, bootstrapped over proteins. This matters: the pooled estimator fails a within-protein
permutation control, returning -0.169 [-0.272, -0.057] on shuffled data, because proteins differ
both in overall indel cost and in class composition. After the change, 200 permutations give a
null of +0.0009 (sd 0.060).

## 3. The deletion penalty, and how it depends on tolerance

Within protein, the mean offset cost(deletion) - cost(insertion) is **0.303 [0.257, 0.351]**
kcal/mol.

It is larger at tolerant positions. The honest way to show this is to bin positions by a *second,
independent* measurement - the glycine insertion at the same position - and then read the offset
between the deletion and the alanine insertion:

| glycine insertion cost (bin mean) | offset, del - Ala ins | 95% CI |
|---|---|---|
| -0.10 | 0.420 | [0.337, 0.507] |
| 0.38 | 0.555 | [0.445, 0.672] |
| 0.89 | 0.515 | [0.418, 0.617] |
| 1.44 | 0.431 | [0.348, 0.519] |
| 1.99 | 0.421 | [0.338, 0.509] |
| 2.59 | 0.331 | [0.245, 0.415] |
| 3.45 | 0.323 | [0.224, 0.422] |
| 5.61 | 0.262 | [0.132, 0.392] |

**Why the independent binning matters.** Binning on the alanine insertion cost itself - the
obvious thing to do - gives a decline from 0.590 all the way to -0.022, suggesting the two edits
become identical at intolerant positions. That is measurement error on the binning axis: positions
with a low *measured* insertion cost are enriched for negative noise, and the independently
measured deletion is then higher. The same artefact deflates the within-protein regression slope
of deletion on insertion cost; instrumenting with the glycine insertion moves it from **0.712
[0.682, 0.740] to 0.870 [0.834, 0.906]**. Sub-proportionality is real but about a third of what
the naive fit implies, and the convergence to equality is not real at all.

Assay ceiling was checked separately: Tsuboyama clip dG to [-1, 5] in their own analyses, and
ddG therefore saturates at roughly wild-type dG + 1. Only 0.9% of insertion and 2.0% of deletion
variants lie within 0.25 kcal/mol of the first-percentile floor here, and dropping every position
where either edit is within 0.5 of it (97% of data retained) leaves the slope and offsets
unchanged.

## 4a. Deletions: all 20 residues, and the interaction has the predicted shape

The deleted residue is the wild-type one, so unlike insertions the deletion side covers every
amino acid: 19,897 single deletions across 387 domains, 20 residue types, 16 of them with at
least 100 deletions in both helices and loops.

**The interaction is real.** An OLS of deletion cost on domain fixed effects, residue identity,
structure class, relative solvent accessibility and residue volume, compared with the same model
plus a residue-by-structure interaction, gives **F = 7.375, p = 1.7e-36** with standard errors
clustered by domain. Permuting the residue identity within domain abolishes it (F = 0.877,
p = 0.68).

**It has the shape helix propensity predicts.** For each residue we take the within-domain mean
deletion cost in helices minus that in loops, and correlate those 19 contrasts with the Pace &
Scholtz helix penalty (alanine 0.0, glycine 1.00, proline 3.16 - higher is a worse helix former):

> **Spearman = -0.603, bootstrap CI over domains [-0.691, -0.446]**

Negative means a *good* helix former costs *more* to remove from a helix, which is what a
propensity account requires. Residue volume gives a weaker and non-significant correlation with
the same contrasts (+0.326, p = 0.30), so this is not size in disguise. The individual contrasts
range from +1.33 (lysine) and +1.19 (alanine) down to +0.05 (valine) and +0.24 (threonine) - the
beta-branched residues, which are the poorest helix formers among the common ones.

The corresponding strand contrasts are positive for 18 of 19 residues (methionine +1.52, alanine
+0.88, down to aspartate +0.06, with asparagine the single negative at -0.41), so deleting from a
strand costs more than deleting from a loop almost regardless of what is deleted - the strand
penalty is much less residue-dependent than the helix penalty.

## 4b. Insertions: two residues, and the sign flips between them

Asymmetry A = cost(insertion) - cost(deletion), within protein then over proteins; negative means
the insertion is the cheaper edit.

| | strand - loop | helix - loop |
|---|---|---|
| **alanine** insertion | **-0.324** [-0.423, -0.227] | -0.082 [-0.150, -0.016] |
| **glycine** insertion | -0.066 [-0.161, 0.026] | **+0.260** [0.203, 0.325] |

The two arms disagree by about 0.26 kcal/mol in the strand-loop contrast and by about 0.34 in the
helix-loop contrast, and change the sign of the latter. Independently, the difference
ddG(Ala) - ddG(Gly) for an inserted residue in this dataset is +0.258 in strands, +0.309 in
helices and -0.017 in loops: the gap between the two arms **is** that preference. In other words,
what looks like a property of beta strands is in large part a property of alanine in beta strands.

**The residue dependence replicates in the designed miniproteins**, with the same direction and a
smaller size:

| designed HEEH_KT | strand - loop | helix - loop |
|---|---|---|
| alanine insertion | -0.192 [-0.275, -0.113] | -0.126 [-0.209, -0.049] |
| glycine insertion | -0.112 [-0.180, -0.045] | **-0.010 [-0.088, 0.065]** |

Switching the inserted residue abolishes the helix contrast there, as it reverses it in natural
domains. The underlying preference also replicates: cost(Gly) - cost(Ala) for an inserted residue
in the designs is +0.208 [0.179, 0.238] in helices, +0.169 [0.116, 0.224] in strands and +0.092
[0.061, 0.122] in loops, against +0.309 / +0.258 / -0.017 in natural domains.

An independent test in Topolska et al.'s yeast data, whose 488 single-residue insertions span all
20 residues, was **inconclusive**: the interaction between the inserted residue's helix propensity
(Pace & Scholtz) and being in a helix is +0.070 [-0.102, 0.242], the right sign but far from
significant (F-test p = 0.67), and a scale-free residue-by-structure interaction gives p = 0.15.
With 112 helix positions spread over 20 residues that test has little power, so it neither
supports nor contradicts.

For the alanine arm the strand excess is robust: it survives restriction to positions whose
wild-type residue is alanine (-0.842 [-1.133, -0.561], so both edits act on the same chemistry),
all three residue-volume tertiles, 13 of 16 wild-type residue types, domain-clustered regression
on volume and burial and pLDDT, the Lehner lab's own structure labels, dropping every position
whose deleted residue is a protease site (-0.276 [-0.381, -0.175]), and - importantly - all six
bins of an independent tolerance measure (-0.507, -0.442, -0.656, -0.318, -0.525, -0.370), so it
is not the tolerance effect of §3 in disguise. In the designed miniproteins the alanine arm gives
-0.192 [-0.275, -0.113] and the glycine arm -0.112 [-0.180, -0.045].

The practical consequence is the title. A statement of the form "strands tolerate indels worse
than helices" has no determinate value: the sign of the helix-loop contrast here flips between two
single-residue insertions that differ by one methyl group.

## 5. The assay is not the explanation

Tsuboyama's dG combines a measured protease-resistance K50 with a sequence model of the unfolded
state, and a beta strand's neighbourhood is hydrophobic, so a biased unfolded-state model could
manufacture a class difference. Strand minus helix, alanine insertion:

| quantity | difference | 95% CI |
|---|---|---|
| dG | +0.557 | [0.330, 0.783] |
| raw K50, trypsin | +0.260 | [0.159, 0.359] |
| raw K50, chymotrypsin | +0.307 | [0.201, 0.414] |
| unfolded-state model, trypsin | +0.003 | [-0.000, 0.005] |
| unfolded-state model, chymotrypsin | +0.003 | [-0.003, 0.009] |

The model contributes about 1%; the raw measurements carry the effect; the two proteases agree
across domains at r = 0.929.

## 6. Three negative results

**(a) Beta strands are not more indel-tolerant than helices.** Paired within protein, natural
domains gave insertions -0.478 [-0.661, -0.288] and deletions -0.445 [-0.669, -0.223] - strands
cheaper, against the received view - and it survived burial, pLDDT, element position, residue
composition and the Lehner lab's labels (-0.547 [-0.731, -0.365] over 134 domains). It does not
replicate: in the designed miniproteins it is -0.024 [-0.111, 0.067], and scaled by each protein's
wild-type stability the designed interval [-0.059, 0.031] excludes the natural estimate of -0.117.
Natural domains restricted to the designs' stability range weaken to -0.278 [-0.555, 0.008]. We
report it because of how it looked before the replication: large, five-way controlled,
assay-verified, and contradicting two published claims.

**(b) The T4 lysozyme accommodation magnitudes do not transfer.** Vetter et al. (1996) report
4.9 kcal/mol for bulging inside a helix against 2.6 for looping out near its end, implying about
2.3 for mid-helix over helix-end. Across 223 domains the paired difference is **0.481
[0.387, 0.575]** for alanine and 0.721 [0.621, 0.823] for glycine, surviving every tertile of
burial and pLDDT. The direction holds; the magnitude is five times smaller. Heinz et al. (1994)
further describe insertions being absorbed by a register shift into the adjacent loop; across 242
domains the correlation between asymmetry and distance to the nearest loop is -0.023
[-0.058, 0.012].

**(c) The beta-bulge mechanism is not established.** Two of three pre-registered marks passed -
the asymmetry is larger in antiparallel than parallel bridges (-0.295 [-0.499, -0.093]) and grows
toward the middle of a strand (-0.304 [-0.504, -0.114]) - but edge strands, which should bulge
most freely, were not reliably different from interior ones (-0.107 [-0.234, 0.022]). Bulges
themselves could not be detected: strict DSSP bridge rules assign no bridge at a bulged residue.
Topolska et al. report bulge accommodation as rare in these same data.

## 7. A note on metrics

Published indel predictors are compared across incompatible statistics: INDELi reports a
per-domain mode Pearson of about 0.62, IFUM a pooled Pearson of 0.80 on the same measurements.
For simple structural features here, pooling across domains **lowers** r (pLDDT 0.383 pooled
against 0.514 per domain) because it injects between-domain variance a within-domain feature
cannot explain, while a model that also predicts domain-level stability gains from that same
variance. The direction of the bias depends on the predictor, so the two numbers are not merely on
different scales - they cannot be ordered.

## 8. Limitations

One assay (protease resistance), two protein sets from one deposit, two inserted residues. An
orthogonal yeast abundance assay on nine domains shows the right ordering (strand -0.224, loop
-0.146, helix +0.003) but failed its pre-registered mark, and only three of those nine domains
contain enough positions of both classes to estimate a contrast, so it is uninformative rather
than contradictory. The mechanism of the residue dependence is untested. All pre-registrations,
including the marks that failed and the two framings withdrawn during the analysis, are deposited
with the code.
