# Indel accommodation: results, 26 Nov 2026

Pre-registrations: `PREREG_indel.md`, `PREREG_bulge.md`, `PREREG_strand.md`.
Literature done first: `LITERATURE_indel.md`.
Every confidence interval is a bootstrap over **domains**, never over positions. Every class
comparison is computed **within domain** first and then averaged, because the pooled form failed
its own permutation control (see §0).

Convention: `cost = dG(wild type) - dG(variant)`, so larger cost = more destabilizing.
`asymmetry A = cost(insertion) - cost(deletion)`; negative A = the insertion is the cheaper edit.

---

## Summary of what survived

| claim | verdict |
|---|---|
| The T4-lysozyme helix accommodation magnitudes generalise | **refuted** - 0.48 kcal/mol where 2.3 was predicted |
| Indel cost depends on distance to the nearest loop (H2) | **refuted** - rho = -0.02 |
| Beta strands tolerate indels better than helices within a protein (H-strand) | **refuted** - fails in designed proteins and in low-stability natural ones |
| Insertions are cheaper than deletions at the same position (0.303 [0.257, 0.351]) | **supported**; the qualitative claim is already published |
| The deletion penalty is larger at tolerant positions | **supported, halved** - it falls 0.52 to 0.26, it does **not** vanish; the vanishing was a regression artefact |
| Insertions are cheapest, relative to deletions, in beta strands | **supported for alanine insertions only**; absent for glycine in natural domains |
| "Indel tolerance by secondary structure" is a well-defined quantity | **refuted** - the helix-loop contrast changes sign between two inserted residues |
| Pooled and per-domain correlations are interchangeable | **refuted** - and the direction of the bias depends on the predictor |
| This project's earlier secondary-structure labels were sound | **refuted** - corrected, Gate 11 recomputed |

---

## 0. A methodological correction that affects everything, including earlier gates

The structure classes used throughout this project came from `ss3()` in
`16_gate11_gly_ala.py`: two very wide Ramachandran boxes that label about 89% of residues as helix
or strand and leave only ~10% as loop. `analysis/17_dssp.py` replaces them with Kabsch-Sander
hydrogen-bond assignment on the AlphaFold models.

Validation, against geometry the algorithm never used (`18_ss_check.py`):

| class | phi median | psi median | CA(i)-CA(i+4) |
|---|---|---|---|
| H | -63.0 | -41.8 | 6.2 A |
| E | -113.7 | +134.1 | 12.2 A |

Those are the textbook values. Against the Lehner lab's own table for the same domains
(`23_ss_crosscheck.py`, 13,682 positions, 239 domains): three-state agreement **0.801**, residue
identity 1.00, rSASA correlation 0.955, element-length correlation 0.941, and **zero** H/E
confusion - the only disagreement is that this assignment ends elements about one residue earlier.

**Consequence for Gate 11** (`19_gate11_dssp.py`, 19,291 positions, 385 domains):

| ddG(Ala) - ddG(Gly) | old labels | hydrogen-bond labels |
|---|---|---|
| helix | +0.201 | **+0.309** |
| strand | +0.102 | **+0.258** |
| loop | -0.082 | **-0.017** |
| helix minus loop | +0.283 [0.244, 0.321] | **+0.326 [0.295, 0.358]** |

Gate 11's pre-registered mark B (>= 0.3) had failed at 0.283; with corrected classes it passes.
**This must always be reported with the reason**: the correction was made as a prerequisite for a
different question, was validated independently, and the marks were not touched - but an analysis
change that turns a failed mark into a passed one deserves the reader's suspicion, and both
numbers belong side by side.

The corrected picture is also cleaner physically: loops are **indifferent** between Ala and Gly
(-0.017, median +0.004); the apparent glycine preference in loops was an artefact of a
contaminated class. New: Ala's advantage is 0.204 [0.162, 0.247] larger in helix middles than at
helix ends.

### The estimator problem

`mean(class A) - mean(class B)` pooled over positions failed its permutation control: on
within-domain-shuffled data it gave **-0.169 [-0.272, -0.057]** instead of zero, because domains
differ both in overall cost and in how many positions of each class they contribute, and a
within-domain shuffle preserves that. Every statistic below is therefore computed inside each
domain first. After the change, 200 independent permutations give a null of **+0.0009 (sd 0.060)**
for insertions and **-0.0037 (sd 0.067)** for deletions, with the observed values outside all 200
draws.

---

## 1. The T4 lysozyme accommodation model does not transfer (H1: failed)

Vetter et al. 1996 report, from 32 insertions in nine helices of T4 lysozyme, that bulging inside
a helix body costs 4.9 kcal/mol against 2.6 for looping-out near a helix end - a predicted
middle-minus-end difference of about 2.3.

Measured here on 223 domains, paired within domain:

> **helix middle minus helix end = 0.481 kcal/mol [0.387, 0.575]** (Ala insertion)
> glycine arm: 0.721 [0.621, 0.823]

- mark H1a (>= 1.0 with CI lower bound > 0.5): **failed**
- mark H1b (CI contains 2.3): **failed**
- mark H1c (same sign in every tertile of rSASA and of pLDDT): **passed**

So the direction is right and robust, and the magnitude is about five times smaller than the
crystallographic estimates. Topolska et al. had reported this difference only qualitatively
("both ends of alpha-helices are indeed slightly more tolerant"); this puts a number on it and
shows that the T4 figures are specific to that protein.

## 2. Loop proximity does not govern indel cost (H2: failed)

Heinz et al. 1994 describe insertions being absorbed by a register shift toward the N terminus
that spills residues into the **preceding loop**. If loop access mattered, the asymmetry inside
helices should depend on distance to the nearest loop. It does not:

> Spearman(distance to nearest loop, A) within helices = **-0.023 [-0.058, 0.012]**

and the asymmetry is the same near a loop (-0.234) and far from one (-0.251). Marks H2a, H2b and
H2c all failed. The local-absorption picture taken from one crystal structure has no detectable
counterpart across 242 domains.

## 3. The surviving result: insertions are cheaper than deletions, most of all in beta strands

Asymmetry A by class, within domain then over domains (Tsuboyama Lib3, Ala insertions):

| class | A (kcal/mol) | 95% CI | domains |
|---|---|---|---|
| beta strand | **-0.587** | [-0.710, -0.468] | 228 |
| alpha helix | -0.266 | [-0.331, -0.200] | 242 |
| loop | -0.243 | [-0.295, -0.192] | 385 |

**The comparison that holds up is strand against loop, paired inside the same protein:**

> **E - L = -0.324 kcal/mol [-0.423, -0.227]** over 227 domains
> H - L = -0.061 [-0.124, 0.003] - helices are not distinguishable from loops
> E - H = -0.033 [-0.171, 0.108] over the 84 domains containing both - not distinguishable

The large strand-versus-helix gap in the unpaired table is mostly **between proteins**: of 387
domains, 159 contain no strand and 145 no helix, and beta-rich domains are globally more
indel-sensitive (Spearman between beta fraction and asymmetry -0.178, p = 0.0005).

### Controls (22_asym_control.py, 26, 29, 31, 32)

- **residue identity**: restricted to positions whose wild-type residue *is* alanine, so the
  deleted and inserted residues are chemically identical, the strand value is **-0.842
  [-1.133, -0.561]** - larger, not smaller
- **volume**: present in all three residue-volume tertiles (-0.868 / -0.451 / -0.683), no trend
- **per residue type**: strand more negative than non-strand in 13 of 16 types, median -0.35
- **regression**: is_E coefficient -0.198 after volume, rSASA and pLDDT, domain-clustered errors
- **burial**: in helices and loops the asymmetry nearly vanishes when exposed (-0.086, -0.153);
  in strands it persists (-0.571)
- **not a scale effect**: the within-domain slope of deletion cost on insertion cost is **0.695
  [0.667, 0.723]**, below 1, and the asymmetry is *largest* where insertion cost is *lowest*
  (-0.758 in the cheapest sixth, -0.265 in the dearest). At matched insertion cost the strand
  excess over helix survives in all six bins, mean -0.384
- **permutation**: 200 draws, null +0.0009 (sd 0.060), observed outside every draw
- **independent structure labels**: holds with the Lehner lab's own assignment
- **not an artefact of the assay's unfolded-state model**: see below

### The assay cannot explain it

Tsuboyama's dG comes from protease resistance plus a sequence-based model of the unfolded state.
A beta strand's neighbourhood is hydrophobic, so a biased unfolded-state model could manufacture a
class difference. It does not (`32_unfolded_check.py`, strand minus helix, Ala insertion):

| quantity | difference | 95% CI |
|---|---|---|
| deltaG | +0.557 | [0.330, 0.783] |
| **raw K50, trypsin** | **+0.260** | [0.159, 0.359] |
| **raw K50, chymotrypsin** | **+0.307** | [0.201, 0.414] |
| unfolded-state model, trypsin | +0.003 | [-0.000, 0.005] |
| unfolded-state model, chymotrypsin | +0.003 | [-0.003, 0.009] |

The unfolded-state model contributes about 1%. The raw measurements carry the effect, and the two
proteases agree across domains at Pearson **0.929**.

### Replication

**Designed miniproteins (passed, and this is the strong one).** Lib3 also contains single indels in
69 designed HEEH_KT miniproteins, 57 with AlphaFold models here, about 131 variants each; neither
Topolska et al. nor Gutierrez & Rocklin appear to have used them. Numbering verified at **1.00**
over 2,412 deletion variants. These proteins are not natural and share no evolutionary history
with the discovery set.

| class | A, designed | 95% CI |
|---|---|---|
| beta strand | **-0.136** | [-0.221, -0.053] |
| alpha helix | -0.055 | [-0.114, 0.000] |
| loop | +0.007 | [-0.046, 0.060] |

> paired **E - L = -0.147 [-0.220, -0.076]** over 53 designs
> (natural domains: -0.327 [-0.422, -0.233])

Same sign, same ordering, in proteins that have no shared ancestry with the first set.

**Yeast abundancePCA (inconclusive, not contradictory).** Topolska et al.'s own orthogonal assay,
9 domains, 491 matched positions. Class means: strand -0.224 [-0.437, -0.087], loop -0.146, helix
+0.003 - the ordering is right. But the pre-registered mark R1 (strand minus non-strand, negative
in >= 7 of 9) **failed** at +0.041, because loops are also negative there and the pooling of helix
with loop dilutes the contrast. R2 (exposed positions only) passed at -0.135 [-0.284, -0.011].
R3 could not be run: there are no insertions whose inserted residue matches the wild type.
Measurement error is comparable between the two edit types (sigma 0.084 against 0.088), so the
failure is not an error artefact.

### CORRECTION, written after the two sections below were first drafted

Two claims made earlier in this file were **withdrawn on controls run afterwards**. Both are kept
with their original numbers, because the reason each failed is more useful than the claim was.

**Withdrawal 1: "the asymmetry vanishes at intolerant positions" is a regression artefact.**
The decline was measured by binning positions on the *same* alanine insertion cost that enters the
offset, so measurement error on the binning axis manufactures a declining line that crosses zero
near the mean. Binning instead on the **glycine** insertion at the same position - an independent
measurement - the offset falls from 0.420/0.555 at the tolerant end only to **0.262
[0.132, 0.392]** at the intolerant end, and **every bin excludes zero**. Likewise the regression
slope of deletion on insertion cost moves from 0.712 [0.682, 0.740] to **0.870 [0.834, 0.906]**
when instrumented with the glycine insertion. Corrected statement: the deletion penalty is about
twice as large at tolerant positions as at intolerant ones, and never disappears.

**Withdrawal 2: the strand excess is specific to the inserted residue.**

| | strand - loop | helix - loop |
|---|---|---|
| alanine insertion | -0.324 [-0.423, -0.227] | -0.082 [-0.150, -0.016] |
| glycine insertion | **-0.066 [-0.161, 0.026]** | **+0.260 [0.203, 0.325]** |

The two arms differ by about 0.26 in the strand contrast and flip the sign of the helix contrast.
That gap **is** Gate 11: ddG(Ala) - ddG(Gly) for an inserted residue is +0.258 in strands, +0.309
in helices, -0.017 in loops. So a large part of what looked like a property of beta strands is a
property of *alanine in* beta strands. In the designed miniproteins the glycine arm does survive
(-0.112 [-0.180, -0.045]) while in natural domains it does not.

What still stands for the alanine arm: the strand excess survives all six bins of the independent
tolerance measure (-0.507, -0.442, -0.656, -0.318, -0.525, -0.370), so it is not the tolerance
effect in disguise.

**The residue dependence itself replicates in the designed miniproteins** (2,348 positions, 57
designs): Ala gives strand-loop -0.192 [-0.275, -0.113] and helix-loop -0.126 [-0.209, -0.049],
Gly gives -0.112 [-0.180, -0.045] and **-0.010 [-0.088, 0.065]** - switching the inserted residue
abolishes the helix contrast there as it reverses it in natural domains. **Gate 11 also replicates
in this independent set**: cost(Gly) - cost(Ala) is +0.208 [0.179, 0.238] in helices, +0.169
[0.116, 0.224] in strands, +0.092 [0.061, 0.122] in loops (natural: +0.309 / +0.258 / -0.017).

An independent attempt in Topolska et al.'s yeast data, whose 488 single-residue insertions span
all 20 residues, is **inconclusive**: helix-propensity-by-helix interaction +0.070 [-0.102, 0.242]
(right sign, F-test p = 0.67) and a scale-free residue-by-structure interaction at p = 0.15. Only
112 helix positions spread over 20 residues - no power.

A correction to the earlier aPCA work in this file: those analyses pooled 1, 2 and 3 residue
indels, which the data distinguishes only through the variant sequence length. Restricted to
single-residue indels (483 matched positions), the asymmetry by class is strand -0.131
[-0.187, -0.081], loop -0.040 [-0.118, 0.040], helix **+0.057 [-0.002, 0.117]**, with paired
E - H = -0.174 [-0.267, -0.085] and E - L = -0.044 [-0.120, 0.036]. The overall deletion penalty
in that assay is +0.043 [-0.025, 0.106] in their scaled abundance units, i.e. essentially absent,
which is a real discrepancy with the +0.303 kcal/mol found in the stability data and is not
resolved here.

### The asymmetry only exists where there is headroom (35_law.py) - SUPERSEDED, see the correction above

Fitting deletion cost against insertion cost **within each protein**: slope **0.692
[0.664, 0.720]**, intercept **0.81 [0.71, 0.92]** kcal/mol. A slope below one means the two edits
converge as positions become less tolerant. By quartile of insertion cost:

| | insertion cost | deletion cost | offset |
|---|---|---|---|
| natural, most tolerant quartile | 0.024 | 0.454 | **+0.554 [0.487, 0.624]** |
| natural, least tolerant quartile | 4.379 | 4.497 | **-0.013 [-0.101, 0.076]** |
| designed, most tolerant quartile | -0.179 | 0.059 | **+0.267 [0.203, 0.334]** |
| designed, least tolerant quartile | 1.156 | 1.091 | **-0.190 [-0.322, -0.062]** |

**Not a measurement ceiling.** Only 0.9% of insertion and 2.0% of deletion variants lie within
0.25 kcal/mol of the 1st-percentile dG floor. Dropping every position where either edit is within
0.5 of that floor keeps 97% of the data and changes nothing: slope 0.698 [0.669, 0.726],
intercept 0.771, tolerant-quartile offset +0.545 [0.481, 0.612], intolerant +0.016 [-0.054,
0.085]. Decisively, in the designed proteins the offset has already **reversed** at a mean
insertion cost of only 1.11 kcal/mol, far from any ceiling.

**The strand effect runs against this pattern, so the two are independent.** Strands are among
the costliest places for an indel, so the headroom law predicts they should show the *least*
asymmetry. They show the most (offset by class: strand +0.596, helix +0.273, loop +0.236), and the
strand excess over helix survives in all six matched insertion-cost bins (mean -0.384).

### Protease-site control

Dropping every position whose deleted wild-type residue is a trypsin (K, R) or chymotrypsin
(F, W, Y, L, M) site - 38% of the data, a precaution prompted by the Rocklin lab excluding
designed EHEE variants as possible protease-site artefacts - leaves the strand-minus-loop
asymmetry at **-0.276 [-0.381, -0.175]** in natural domains and **-0.105 [-0.189, -0.022]** in
designed ones. About a third smaller, still clearly non-zero.

### Mechanism (PREREG_bulge.md): one of the two load-bearing marks passed

H-bulge was that a beta bulge gives an inserted residue an accommodation route that a deletion
does not have.

| mark | result |
|---|---|
| M1 antiparallel more asymmetric than parallel | **passed**: -0.295 [-0.499, -0.093]; anti -0.633, parallel -0.438 |
| M2 edge strands more asymmetric than interior | **failed**: -0.107 [-0.234, 0.022] (edge -0.635, interior -0.461) |
| M3 asymmetry grows toward the strand middle | **passed**: -0.304 [-0.504, -0.114] |

Pre-registered verdict for one of two load-bearing marks: **partial**. The decomposition is
consistent with the hypothesis in direction - relative to each domain's own loops, the *insertion*
cost in strands is only +0.271 [0.175, 0.368] against +0.692 [0.612, 0.773] in helices, so the
asymmetry comes from insertions being cheap rather than deletions being dear - but M2 failed and
the bulge itself could not be detected: strict DSSP bridge rules assign no bridge at a bulged
residue, so 0 bulges were found in 862 models and the "already bulged" test is not possible with
this tooling. **H-bulge is not established.**

## 4. Strands are NOT more indel-tolerant than helices (H-strand: refuted)

In natural domains the within-protein paired difference looked strong and passed all five
controls - burial, pLDDT, element position, residue composition, and the Lehner lab's labels
(-0.547 [-0.731, -0.365] over 134 domains with their assignment):

> natural, insertion: -0.478 [-0.661, -0.288]; deletion: -0.445 [-0.669, -0.223]

It does not replicate.

| | natural | designed |
|---|---|---|
| strand minus helix, insertion | -0.478 [-0.663, -0.290] | **-0.024 [-0.111, 0.067]** |
| same, scaled by the protein's wild-type dG | -0.117 [-0.185, -0.047] | **-0.015 [-0.059, 0.031]** |

The scaled designed interval excludes the natural estimate, so this is not merely a dynamic-range
problem. And within natural domains restricted to the designs' stability range (266 domains), the
effect weakens to -0.278 [-0.555, 0.008] for insertions and disappears for deletions
(-0.148 [-0.446, 0.130]).

**Conclusion: the strand-versus-helix difference is tied to stable natural domains and is not a
general property of beta strands.** It is withdrawn. Note what this would have looked like without
the designed-protein check: a large, five-way-controlled, assay-verified effect that contradicts
two published claims.

## 5. Metric audit (B1: failed, with the sign reversed)

INDELi reports a per-domain mode r of 0.62; IFUM reports a pooled Pearson of 0.80 on the same
underlying data. For simple structural features here:

| predictor, Ala insertions | pooled abs r | per-domain median abs r |
|---|---|---|
| pLDDT | 0.383 | **0.514** |
| rSASA | 0.233 | 0.285 |
| distance to terminus | 0.207 | 0.302 |

Mark B1 predicted pooling would **inflate** r by >= 0.10; it **deflates** it by about 0.10. The
reason: pooling adds between-domain variance in cost, which a within-domain feature cannot
explain, so it dilutes r - while a model that also predicts domain-level stability (IFUM does)
gains from exactly that variance. **The direction of the bias depends on the predictor**, which
makes a pooled number and a per-domain number not merely differently scaled but not orderable at
all. The literature is currently comparing them.

---

## Files

`17_dssp.py` structure, `18_ss_check.py` geometric validation, `19_gate11_dssp.py` Gate 11 redo,
`21_indel_marks2.py` marks B/H1/H2/H3, `22_asym_control.py` identity controls,
`23_ss_crosscheck.py` external label check, `25_bridges.py` bridge partners, `26_mechanism.py`
M1-M3, `27_within_domain.py` within vs between, `28_decompose.py` components,
`29_multiplicative.py` scale control, `30_strand_controls.py` C1-C5 and aPCA,
`31_adjudicate.py` permutation stability and shared domains, `32_unfolded_check.py` assay check,
`33_designed.py` designed-protein replication, `34_scale.py` common-scale comparison.

---

## 6. The deletion side: all 20 residues (PREREG_deletion_residue.md, 39_deletion_residue.py)

The residue-dependence claim above rested on one contrast, alanine against glycine, because those
are the only two residues inserted anywhere in the Tsuboyama libraries (50,049 Ala, 50,139 Gly,
no others). The deletion side removes the wild-type residue, so all 20 types are present and the
same question can be asked with far more power: 19,897 deletions, 387 domains, 16 residues with
at least 100 deletions in both helices and loops.

| mark | result |
|---|---|
| **D-1** residue-by-structure interaction, after domain, burial and volume | **passed**: F = 7.375, p = 1.7e-36, domain-clustered |
| kill 1, permute residue identity within domain | **vanishes**: F = 0.877, p = 0.68 |
| **D-2** helix contrast scales with Pace helix propensity | **passed**: Spearman **-0.603**, bootstrap CI over domains **[-0.691, -0.446]** |
| **D-4** propensity beats volume | **passed**: volume gives +0.326, p = 0.30 |
| D-3 within burial tertiles | not evaluable - too few residues per stratum |
| D-5 designed replication | underpowered: 8 residues, Spearman -0.204, p = 0.63 |

Note on D-2: the first pass computed a Spearman p-value over the 12 residue points and got
-0.539, p = 0.07, which fails. The **pre-registered** test was a bootstrap over domains, and that
gives -0.603 [-0.691, -0.446] over 19 residues. The pre-registered version is the one reported;
the weaker first pass is recorded here so the difference is visible.

Helix-minus-loop deletion cost by residue (positive = costlier to delete from a helix):
K +1.33, Q +1.28, R +1.25, A +1.19, D +1.13, I +1.05, E +0.98, S +0.57, L +0.41, T +0.24,
N +0.23, V +0.05. Strand-minus-loop is positive for 18 of 19 residues (M +1.52 down to D +0.06,
N -0.41), so the strand penalty is far less residue-dependent than the helix penalty.

**What this does for the central claim.** The insertion arm shows that switching one inserted
residue reverses the helix contrast; the deletion arm shows, with 20 residues and ten times the
statistical power, that the structure dependence is an ordered function of the residue's helix
propensity. Together they say the secondary-structure dependence of indel cost is a property of
the residue involved rather than of the indel.

---

## Addendum, 2026-11-26: self-verification before external review

Two things were checked that had been left open, and one of them changed the manuscript.

### Proline (closed, no impact)

The 19 helix propensities are the values listed in the abstract of Pace & Scholtz (1998),
each checked against that source. The full text sits behind an anti-automation
interstitial on PMC and was not retrieved, so a proline value could not be verified; a
figure of 3.16 kcal/mol circulates in secondary sources.

Rather than chase the number, `analysis/48_proline_sensitivity.py` tested whether it
matters. It does not: the correlation is rank-based and proline ranks last on every
account, so giving it any value from 1.0 to 5.0 moves the correlation only from -0.64 to
between -0.59 and -0.58, and every value at or above 1.5 returns an identical figure.

An incidental finding worth keeping: proline is the one residue that does not follow the
trend. Its measured helix-minus-loop contrast is +0.61 kcal/mol over 68 proteins, placing
it 10th of 19 from the most negative end, where the propensity account predicts the
extreme alongside glycine (-0.05). The manuscript now reports this instead of hiding it
inside an exclusion. The earlier sentence "the published scale excludes it" was removed:
the abstract's phrasing ("excluding proline, glycine has the lowest") is ambiguous and
that claim could not be supported.

### Secondary structure, third independent check (changed the manuscript)

The assignment had been validated twice, but both checks were computed from the same
AlphaFold coordinates. `analysis/49_pdb_crosscheck.py` removes that shared dependency: each
domain is named by a PDB identifier, so each has a deposited experimental structure whose
HELIX and SHEET records were assigned by the depositors from different coordinates.

Over 410 domains and 22,738 residues (447 entries, mostly NMR, carry no such records):
three-state agreement **0.802**, helix/strand confusions **7 in total (0.03%)**. The 0.802
is independently almost identical to the 0.801 obtained earlier against a different
external table.

The disagreement is one-sided: **3,866** positions where we say loop and the deposited
record says helix or strand, against **619** the other way. Per-class agreement with
experiment is 0.96 for helix, 0.90 for strand, **0.67 for loop**.

**The prediction that followed was wrong, and this is the part that matters.**
`analysis/50_consensus_ss.py` declared in advance that cleaning the loop class should
*increase* the contrast, making our estimate conservative. Half of that held and half did
not:

| variant | mean contrast | Spearman vs propensity |
|---|---|---|
| our labels, 394 covered domains | 0.76 | -0.54 [-0.67, -0.44] |
| consensus positions only | 0.84 | **-0.35** [-0.55, -0.23] |
| deposited labels substituted | 0.72 | **-0.25** [-0.41, -0.15] |

The contrast grew as predicted. The **correlation with helix propensity, which is the
manuscript's central claim, weakened substantially** — from -0.54 to -0.25 under the
harshest variant. Every interval still excludes zero, so the ordering is real under all
three assignments, but its strength is not independent of how structure is defined.

Part of this is expected attenuation: depositor annotations are coarser, come from a
different structure, and noise in a categorical predictor pulls a residue-level
correlation toward zero. That explanation cannot be used to dismiss the result, because it
cannot be separated from the possibility that our assignment flatters the effect.

The manuscript now reports the full range in the abstract, in a dedicated Results
paragraph and in the Limitations, stating that the rank correlation lies between -0.25 and
-0.64 depending on the structure annotation and is negative throughout. It no longer
quotes -0.64 unqualified.
