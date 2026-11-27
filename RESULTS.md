# insertion-lab: is predicting permissive insertion sites actually solved?

Question: when you want to put something into the middle of a protein (a domain, or a short tag
such as the split-GFP GFP11 fragment), which positions tolerate it? The user has done this by hand
with GFP11 and wanted to know whether the computational version is an open problem.

## State of the field (searched 25 Nov 2026)

- A 2026 preprint on CRISPR-guided PAM scanning states that permissive insertion sites "are often
  sparse and hard to predict from sequence, structure, or computational models alone".
- Two published studies disagree on the basics. Mathony et al. 2023 (Adv. Sci. 10:e2303496) report
  sequence conservation statistics as the most powerful predictors of domain insertion success;
  a bioRxiv study of ion channels reports that permissibility correlates with dynamic features and
  "does not correlate well with static and conservation-related protein properties".
- One model exists: **ProDomino** (Wolf, Mathony, Niopek; bioRxiv 2024.12.04.626757, now in PMC).
  ESM-2 3B embeddings plus a two-layer perceptron, sequence only, no structure. Trained on 174,872
  natural sequences carrying intradomain insertions (InterPro/CATH-Gene3D intersected with
  UniRef50), where the positive label is the natural insertion site plus its two neighbours and the
  negative is a randomly drawn other position. Weights are released. Licence CC BY-NC 4.0, so it
  cannot be used commercially.

**What is missing, and is the reason for this project.**

1. **No published comparison against simple baselines.** ProDomino reports AUROC 0.84 on AraC and
   0.71 on Cas9 transposon data. It never reports what conservation alone, solvent accessibility
   alone, or secondary structure alone achieve on the same data. Those features were computed, but
   only to show the model's predictions do not correlate with them.
2. **The AraC number is not a clean held-out test.** The authors chose the 1,500-step checkpoint
   after watching performance drop at 3,000 steps on that same AraC benchmark, and write that "it
   was critical to benchmark our model on real-world experimental data early on to refine the
   training regime."
3. **Three of the four available proteins were never used.** The same lab's DI_screen contains
   measured data for AraC, Flp, SigF and TVMV. Only AraC appears in the ProDomino paper.
4. **Short tags are out of the training distribution.** The inserts it learned from have a median
   length of about 120 amino acids. GFP11 is 16. Short peptide tags are never mentioned.

## Gate 0, the data census: passed (25 Nov 2026)

`analysis/00_census.py`, output `results/census.json`. Source: github.com/Niopek-Lab/DI_screen
(per-read tables from the DIP-seq pipeline; no licence is stated in that repository, which matters
if anything here is ever published).

Sites counted are distinct in-frame insertion positions with at least 10 reads in both the
unsorted input library and a sorted gate.

| Protein | Inserted domains with data | Best usable sites |
|---|---|---|
| AraC | ERD, LOV, PDZ, UniRapR, eYFP | 270 |
| Flp | PDZ | 269 |
| SigF | PDZ | 225 |
| TVMV | PDZ | 208 |
| **Total** | 8 protein-domain pairs | **972** |

Also present in the repository and needed for the baselines: reference multiple sequence alignments
for all four proteins (`MSA/*.afa`) and AlphaFold models (`*_AF.pdb`), plus some experimentally
relevant complete structures. T7 has a structure and a sequence but no insertion CSVs, so it is out.

This is far too little to train a general model, and enough to benchmark one.

## Gate 1: pre-registered benchmark (written before any score was computed)

**Question.** Does ProDomino predict measured domain-insertion tolerance better than simple
structural and evolutionary baselines, on proteins its authors did not tune it on?

**Label.** Per site and per inserted domain, enrichment = log2 of (share of reads at that site in
the sorted gate) over (share in the input library). A site is positive if its enrichment is above
zero, which is the definition the ProDomino authors used for their AraC benchmark ("positive
enrichment scores for at least one of the five insert domains"). Sites with fewer than 10 reads on
either side are dropped. For AraC the five domains are pooled by the same "any domain" rule; the
other three proteins have PDZ only.

**Predictors compared**, all scored per position:

- **ProDomino**, released checkpoint, default settings.
- **Conservation**, Shannon entropy of the aligned column from the provided MSA.
- **Solvent accessibility**, relative SASA on the AlphaFold model (Shrake-Rupley).
- **Secondary structure**, loop versus helix or strand (DSSP on the AlphaFold model).
- **pLDDT** of the AlphaFold model, as a cheap proxy for local flexibility.
- **Distance to the nearest terminus**, in residues.
- **Combined baseline**, logistic regression on the five features above, fitted leave-one-protein-out
  so it never sees the protein it is scored on.

**Primary metric.** AUROC per protein.

**Primary pass marks, decided now.** AraC is excluded from the primary test because the ProDomino
authors used it to pick their checkpoint. The primary test is on **Flp, SigF and TVMV**. ProDomino
is judged better than the baselines only if both hold:

1. it beats the combined baseline's AUROC on at least 2 of those 3 proteins;
2. its AUROC pooled over those 3 proteins is higher than the combined baseline's.

**Reported without a pass mark.** AraC, including whether AUROC 0.84 reproduces; each single
feature's AUROC per protein; whether tolerance is domain-specific or site-intrinsic, using AraC's
five inserted domains (correlation of enrichment between domain pairs at the same sites); and
calibration of ProDomino's low scores, since its authors note "low scores are not always accurate".

**Limits fixed in advance.** Four proteins from one laboratory and one assay family, so a result
here is about this benchmark, not about every protein. Three of the four have only PDZ inserted,
so "insertion tolerance" is partly "PDZ tolerance". Binarising a continuous enrichment at zero is
the authors' choice, inherited for comparability, and is also reported as a continuous Spearman
correlation. The read-depth cut of 10 is a choice, and the whole analysis is repeated at 30 as a
robustness check. Whether these four proteins' relatives sat in ProDomino's training set is not
known; its CATH-superfamily split makes overlap less likely but does not exclude it, and this is
recorded as an unresolved confound rather than a settled fact.

**What each outcome would mean.** If ProDomino clearly beats the baselines on untouched proteins,
the problem is in better shape than the field's own statements suggest, and the useful open question
narrows to short tags. If it does not, then the one published model is not demonstrably better than
conservation plus solvent accessibility, and the problem is open in the plainest sense.

### Two deviations from the pre-registration, declared before any score was computed

Both were found by reading the authors' own notebook (`DI_screen/analysis/1-sequencing_processing.ipynb`)
after writing the plan above.

1. **The label follows their pipeline exactly, not mine.** Residue index is `insertion_site / 3 - 1`
   (they subtract one for the linker serine after the start methionine); reads are kept when
   `in_frame_insertion` and `forward_insertion` are true and `linker_seq` is empty; counts per site
   are turned into a share of all reads in that sample, the enrichment is the ratio of the sorted
   share to the input share, and the score is its log2, with sites absent from the input set to
   missing and `-inf` set to -10. `12` and `22` are replicates, so the score is their mean, which is
   what the authors use throughout. A site is positive if that mean is above zero. **The authors
   apply no read-depth cut**, so the primary analysis applies none either, and the depth cuts of 10
   and 30 from the plan above become robustness checks.
2. **Secondary structure comes from backbone dihedrals, not DSSP**, because no DSSP binary is
   available here. Helix, strand and loop are assigned from phi and psi. Relative solvent
   accessibility is computed with an own Shrake-Rupley implementation rather than Biopython, which
   is not installed.

## Baseline results (25 Nov 2026), before ProDomino was run

`analysis/01_features.py` and `analysis/02_baselines.py`; outputs `results/features.csv`,
`results/baselines.json`, `results/baselines_sites.csv`. 1,164 sites with a score: AraC 290 (79
positive, pooled over its five inserted domains), Flp 414 (30), SigF 244 (48), TVMV 216 (30).

**Single features, AUROC, with the orientation fixed in advance** (more variable, more exposed, in
a loop, lower pLDDT, nearer a terminus are all expected to be more tolerant; below 0.5 therefore
means the expectation was wrong, not that the feature is useless):

| Feature | AraC | Flp | SigF | TVMV | mean |
|---|---|---|---|---|---|
| conservation (Shannon entropy) | **0.828** | 0.605 | 0.620 | 0.509 | 0.640 |
| relative solvent accessibility | 0.673 | 0.632 | 0.481 | 0.745 | 0.633 |
| pLDDT (lower is better) | 0.779 | 0.749 | 0.460 | 0.676 | 0.666 |
| distance to nearest terminus (nearer is better) | 0.566 | **0.882** | 0.351 | 0.597 | 0.599 |
| in a loop | 0.516 | 0.538 | 0.493 | 0.499 | 0.512 |

**Combined baseline, leave one protein out**: AraC 0.704, Flp 0.439, SigF 0.476, TVMV 0.583;
pooled over the three untouched proteins **0.544**, pooled over all four 0.606.

**Three things follow, none of which needed a model.**

1. **The published disagreement is explained.** Mathony et al. report conservation as the strongest
   predictor and the ion-channel study reports that conservation does not correlate. On AraC
   conservation gives 0.828; on TVMV it gives 0.509, which is chance. Both papers are right about
   the protein each looked at.
2. **No feature keeps its direction across proteins.** Distance to the terminus gives 0.882 on Flp
   and 0.351 on SigF, i.e. it reverses. The same happens, less sharply, for pLDDT and accessibility.
3. **So a model fitted on some proteins mispredicts the next one.** The combined baseline reaches
   0.544 on held-out proteins and falls *below* chance on two of them, because it learns the sign
   of each feature from proteins where that sign is different. Any claim about insertion tolerance
   measured within one protein should be assumed not to transfer until shown otherwise.

**Not yet done, and the reason this is not a verdict:** ProDomino has not been run. Its published
AraC AUROC is 0.84 and conservation entropy alone reaches 0.828 here, but my label may differ in
detail from theirs, so no comparison is claimed until the model is scored on these same labels.

## Gate 1 result (25 Nov 2026): the pre-registered test fails, and the reason is informative

`analysis/03_prodomino.py` (released checkpoint, ESM-2 3B layer 36, no tuning) and
`analysis/04_compare.py`. Outputs `results/prodomino.csv`, `results/compare.json`.

**Indexing was checked first.** Scoring at residue offsets -2 to +2 gives the best or near-best
AUROC at offset 0 for three of the four proteins, so the residue convention matches and no
off-by-one artefact is driving anything below.

| Protein | ProDomino | Combined baseline (leave-one-protein-out) | Best single feature, chosen after the fact |
|---|---|---|---|
| AraC (used by its authors to pick their checkpoint) | **0.877** | 0.704 | conservation 0.828 |
| Flp | 0.650 | 0.439 | distance to terminus 0.882 |
| SigF | 0.673 | 0.476 | conservation 0.620 |
| TVMV | 0.646 | 0.583 | accessibility 0.745 |
| **pooled over the three untouched proteins** | **0.495** | 0.544 | - |

**Pass marks: one met, one missed, so Gate 1 fails.** ProDomino beats the combined baseline on
3 of 3 untouched proteins (mark 1 met), but its pooled AUROC over those three is 0.495 against the
baseline's 0.544 (mark 2 missed). Both were required.

**The pooled number is not a ranking failure, it is a calibration failure, and this matters in
practice.** ProDomino beats the baseline inside every single protein yet loses when they are pooled.
Its mean score per protein is 0.478 for Flp, 0.142 for TVMV and 0.074 for SigF, while the share of
tolerant sites runs the other way: 7.2% for Flp, 13.9% for TVMV, 19.7% for SigF. It gives its
highest scores to the protein with the fewest tolerant sites. So the scores rank positions well
*within* a protein and carry almost no meaning *between* proteins. A user holding a new protein and
a score of 0.5 cannot tell whether that is a good site, which is exactly the situation the tool is
meant to serve.

**The AraC number reproduces, and that is the smaller story.** Published 0.84, measured here 0.877.
But Shannon entropy of one alignment column reaches 0.828 on the same protein and labels. On the
protein its authors reported, a 3-billion-parameter model is about four points of AUROC ahead of
counting letters in an alignment.

**On unseen proteins it is much weaker than its headline suggests**: 0.646 to 0.673, against 0.877
on AraC.

**An important caveat against over-reading the last column.** "Best single feature" is picked after
seeing the answer, and it is a different feature for each protein, so it is an optimistic bound and
not a method anyone could use in advance. The fair comparison for a fixed model is the combined
leave-one-protein-out baseline, and ProDomino beats that everywhere. The honest summary is that
**ProDomino is better than any transferable baseline built here, considerably weaker on unseen
proteins than its published figure, and its scores are not comparable across proteins.**

**What is now established in this project**, with the data and code in this folder:

1. The literature disagreement about conservation is a protein-to-protein effect, not a real
   contradiction (0.828 on AraC, 0.509 on TVMV).
2. Simple structural and evolutionary features do not keep their direction across proteins, so a
   model fitted on some proteins can fall below chance on the next (0.439 and 0.476 here).
3. The one published model generalises to unseen proteins at AUROC about 0.65, not 0.84, and its
   scores are not calibrated between proteins.

None of these three is stated in the papers I read, and each follows from public data and a few
hours of compute.

## Gate 2: does insertion tolerance depend on the size of the insert? (written 25 Nov 2026, before computing)

Why this is the question. The practical problem the user cares about is placing a 16-residue GFP11
tag. ProDomino learned from inserts whose median length is about 120 residues. Topolska et al.
measure 1-3 residue insertions and report that relative solvent accessibility predicts substitution
tolerance well (mode r 0.644) but insertion tolerance poorly (0.168), while in the domain-insertion
data analysed above that same feature reaches AUROC 0.745 on TVMV. Those numbers hint that the
rules differ with insert size, but no one has compared sizes at matched positions.

**The data that allows it.** MaveDB `urn:mavedb:00000660-a-1`, "DIMPLE mKir2.1 surface expression"
(Macdonald et al., Genome Biology 2023): 11,055 scored variants in mouse Kir2.1, of which 1,280 are
insertions placed at 434 positions in three sizes at the same sites - **Gly (1 residue, 427
positions), GlySer (2, 426) and GlySerGly (3, 427)** - plus 1,270 deletions and 8,504 substitutions.
Downloaded to `data/mavedb/kir21_scores.csv`. The phenotype is surface expression measured by
antibody staining and sorting.

The matching domain-insertion data for the same protein (Coyote-Maestas et al., Nat Commun 2019)
is **not public**: its statement is "available from the corresponding author upon reasonable
request", with only raw reads in SRA under PRJNA506141. So a true 1-vs-120-residue comparison
inside one protein is not possible from public data, and is recorded here as the experiment that
would settle it.

**Pre-registered questions and marks.**

- **Q1, size dependence within a protein.** Spearman correlation between the 1-residue and the
  3-residue insertion score at the same position in Kir2.1. Reported with a 95% bootstrap interval.
  Marks fixed now, for interpretation rather than pass/fail: above 0.8 means size barely matters
  over this range; below 0.5 means the ordering of good sites changes with size, which would mean
  "insertion tolerance" is not one property.
- **Q2, does a domain-insertion model transfer to short inserts?** ProDomino scores for Kir2.1
  against each of the three insertion sizes, as Spearman and as AUROC using the dataset's own
  threshold for a strongly reduced score. **Pass mark: ProDomino reaches AUROC of at least 0.65 for
  the 1-residue insertions**, which is the level it reached on the held-out domain-insertion
  proteins in Gate 1. Failing that means a model trained on domain insertions does not describe
  short-tag insertions, which is the case that matters for GFP11.
- **Q3, baselines on the same labels.** Relative solvent accessibility, pLDDT, secondary structure
  and distance to the nearest terminus, taken from the AlphaFold model of mouse Kir2.1, scored the
  same way, with the same orientations as before. Conservation is omitted unless an alignment is
  built, and its absence is noted rather than worked around.

**Reported without a mark.** Whether the sites that tolerate 3 residues are a subset of those
tolerating 1; whether the size effect is concentrated in particular secondary structures; and the
correlation between insertion and deletion tolerance at the same position, which the Kir2.1 data
also supports.

**Limits fixed in advance.** One protein, one assay (surface expression, which confounds folding,
trafficking and stability), inserts of one chemical character (glycine and serine, deliberately
flexible, so these are the easiest possible inserts), and a maximum size of 3 residues, still far
from GFP11's 16. A negative result at 3 residues is therefore strong evidence about short inserts;
a positive result would not prove that 16 residues behave the same way.

### Clarification to Gate 2, declared before any prediction was scored

The pre-registration said AUROC would use "the dataset's own threshold". The deposit has no such
threshold (it contains a single synonymous variant), so the rule is fixed here instead:
**Spearman correlation is the primary, threshold-free measure**, and AUROC is secondary, using as
cutoff the midpoint between the two components of a two-component Gaussian mixture fitted to all
insertion scores. The mixture is fitted to the measured scores only and never sees any prediction.

Also recorded now: the screened construct is not wild-type Kir2.1. Aligning the sequence
reconstructed from the variant names to the AlphaFold model of mouse Kir2.1 (UniProt P35561, model
v6) shows a **FLAG tag, DYKDDDDK, inserted after position 115**, in the extracellular loop used for
surface staining (430 of 436 residues match, the 6 mismatches being the tag itself). Numbering maps
to the structure as: position p for p <= 115, p - 8 for p >= 124, with 116-123 the tag and no
structural counterpart. ProDomino is run on the 436-residue construct, so its positions line up
with the measurements directly.

One observation from the measured scores alone, before any model was involved: median insertion
score is -0.08 for one residue, -0.05 for two and -0.42 for three, against -1.06 for deletions.
Deletions are the most disruptive, which reproduces the DIMPLE paper's own conclusion, and three
residues are already measurably worse than one.

## Gate 2 result (25 Nov 2026): the model transfers, and a free baseline beats it

`analysis/06_kir_size.py`, output `results/kir_size.json`. 434 positions in Kir2.1, scores for
insertions of 1, 2 and 3 residues and for 1-residue deletions. Mixture cutoff -1.324 (component
means -3.384 and 0.736); by it, 34% of 1-residue and 42% of 3-residue insertions are disruptive.

**Q1, size dependence.** Spearman between sizes at the same position: 1 vs 2 residues 0.810
(0.772-0.839), 2 vs 3 0.842 (0.812-0.863), **1 vs 3 0.786 (0.746-0.821)**. Insertions against
deletions at the same position: 0.660 and 0.693, so insertion and deletion tolerance are related
but clearly distinct properties. Median scores: -0.08, -0.05, -0.42 for 1, 2, 3 residues against
-0.56 for deletions, reproducing the DIMPLE paper's own finding that deletions hurt most.

By the bands fixed in advance, 0.786 falls just under "size barely matters". **But the aggregate
number hides where it matters.** Of the worst quarter of positions by their own threshold, only
half are shared between 1 and 3 residues (Jaccard 0.507, 72 of 142); at the worst half the overlap
is 0.733. So size changes little in the middle of the distribution and a great deal in the tail,
which is exactly the part anyone choosing a tag site cares about.

**Q2, does a domain-insertion model transfer to short inserts? Yes, the pass mark is met.**
ProDomino reaches AUROC 0.747 for 1-residue insertions (mark was 0.65), 0.766 for 2 and 3, and
0.725 for deletions, with Spearman 0.456 to 0.509. A model trained only on inserts of about 120
residues does carry real signal down to a single residue.

**Q3, and this is the finding: AlphaFold's own confidence score beats it on every label.**

| Predictor | ins 1 | ins 2 | ins 3 | del 1 |
|---|---|---|---|---|
| **pLDDT (lower is more tolerant)** | **0.850** | **0.875** | **0.894** | **0.816** |
| relative solvent accessibility | 0.754 | 0.762 | 0.785 | 0.765 |
| ProDomino | 0.747 | 0.766 | 0.766 | 0.725 |
| distance to nearest terminus | 0.687 | 0.697 | 0.700 | 0.692 |
| in a loop | 0.562 | 0.563 | 0.575 | 0.555 |

pLDDT is a number AlphaFold already prints. It needs no training, no GPU and no model, and it is
10 to 13 AUROC points ahead of the 3-billion-parameter model on this protein, for every insert
size and for deletions.

**Caveats, stated plainly.** This is one protein, and a membrane protein whose assay is surface
expression, so low-confidence regions are disordered loops that both tolerate insertions and are
easy for AlphaFold to flag. The effect may be weaker in compact soluble domains: Topolska et al.
report accessibility predicting insertion tolerance poorly (mode r 0.168) in small domains, while
here it reaches Spearman 0.46 to 0.53. Whether pLDDT keeps this lead across proteins is the
obvious next test and is not answered here. ProDomino was run on the screened construct including
its internal FLAG tag, so its positions match the measurements exactly.

## Gate 3: does the pLDDT result hold on proteins that did not generate it? (written 25 Nov 2026, before computing)

Gate 2 found that AlphaFold's pLDDT beats ProDomino for short insertions in Kir2.1 by 10 to 13
AUROC points. **Kir2.1 generated that claim, so it cannot test it.** Kir2.1 is also a membrane
protein scored by surface expression, the case most favourable to a disorder-based predictor.

**Held-out proteins, both soluble, both from MaveDB, neither seen before this test:**

- **TEM-1 beta-lactamase** (`urn:mavedb:00000100-a-1`), 4,780 single-codon insertions at 276
  positions, so many different inserted residues per position; phenotype is enzyme fitness.
- **Human DHFR** (`urn:mavedb:00001219-a-1..4`), 186 single-residue insertion positions measured in
  four conditions: 30 and 37 degrees, each with and without methotrexate; phenotype is folding,
  via a yeast folding sensor.

**Fixed now.** For TEM-1 a position's score is the **median over all residues inserted there**, and
the identity question is reported separately without a mark. For DHFR the primary score is the
**mean of the four conditions**, with each condition also reported alone. Thresholds come from the
same two-component Gaussian mixture rule as Gate 2, fitted per dataset to measured scores only.
Structures are the AlphaFold models of UniProt P62593 (TEM-1) and P00374 (DHFR), with numbering
mapped by matching the residue identities in the variant names, as was done for Kir2.1. Feature
orientations are unchanged from earlier gates.

**Primary mark.** The claim is supported only if **pLDDT's AUROC is at least ProDomino's on both
TEM-1 and DHFR**. Winning on one is reported as partial and not as support. Losing on both retires
the claim, and Gate 2's result is then recorded as specific to Kir2.1 or to membrane proteins.

**Reported without a mark.** Relative solvent accessibility, distance to terminus and loop status;
Spearman as well as AUROC; whether the identity of the inserted residue matters in TEM-1; and
whether TEM-1's signal peptide (roughly residues 1-23) behaves differently, reported with and
without it.

**Limits fixed in advance.** Three proteins and three different assays (surface expression, enzyme
fitness, folding sensor), so a consistent result across them is evidence about insertion tolerance
in general, while an inconsistent one may be about the assays rather than the proteins. All inserts
here are one residue, so nothing in this gate speaks directly about a 16-residue tag.

## Gate 4: is insertion tolerance a property of segments rather than residues? (written 25 Nov 2026, before computing)

Every predictor in this field, ProDomino and INDELi included, and every baseline used above, scores
a single position. Physically that cannot be the whole story: an inserted stretch has to go
somewhere, so what should matter is not how flexible one residue is but **how much room the
flexible segment containing it has**.

**Hypothesis H4.** Insertion tolerance is a property of the contiguous flexible segment a position
belongs to, and the segment length that matters grows with the size of the insert.

**The feature.** For each position, `run_len` = the length of the contiguous run of positions whose
pLDDT is below a threshold that contains it (0 if the position itself is above the threshold), and
`run_slack` = the number of residues in that run on either side of the position. The pLDDT
threshold is set per protein at its own median, fixed now, so no tuning is possible; the analysis
is repeated at the 40th and 60th percentile as a robustness check. A loop-based version using the
dihedral secondary structure is computed the same way and reported alongside.

**Predictions, each with a mark fixed now.**

- **P1.** `run_len` has a higher AUROC than per-residue pLDDT for insertion tolerance. Mark: true on
  at least 3 of the 4 short-insertion datasets available (Kir2.1 at 1, 2 and 3 residues, TEM-1,
  DHFR; Kir2.1's three sizes count as one dataset at its largest size for this mark).
- **P2, the one that matters.** The advantage of `run_len` over per-residue pLDDT **grows with
  insert size**. Mark: within Kir2.1, (AUROC of run_len minus AUROC of pLDDT) is larger for
  3-residue insertions than for 1-residue insertions, and larger again for the domain insertions in
  the four DI_screen proteins than for any single-residue dataset. A flat or shrinking advantage
  refutes H4 even if P1 passes, because P1 alone could just mean that smoothing helps.
- **P3.** A position's size sensitivity, defined as its 1-residue score minus its 3-residue score in
  Kir2.1, is predicted by `run_slack`: Spearman below -0.2 with a bootstrap interval excluding
  zero, i.e. more slack means less loss when the insert grows.

**What refutes it.** If `run_len` does no better than per-residue pLDDT (P1 fails), or if its
advantage does not grow with insert size (P2 fails), the segment idea is wrong as stated and will
be reported as such. P2 is the real test; P1 passing alone is a smoothing artefact.

**Limits fixed in advance.** pLDDT is a confidence score, not a measurement of flexibility, and
AlphaFold was trained on crystal structures, which already encode which regions are mobile; a
positive result says the segment picture is useful, not that pLDDT measures physics. Run length and
solvent accessibility are correlated, so the combined model is reported too. The domain-insertion
proteins and the short-insertion proteins are different proteins measured in different assays, so
the cross-size part of P2 compares across systems and is the weakest link in the design.

## Gate 5: does the best residue to insert depend on the site? (written 25 Nov 2026, before computing)

TEM-1 (`urn:mavedb:00000100-a-1`) carries 4,780 single-codon insertions over 276 positions, which
is roughly 17 different inserted residues at each position. Published analyses pool them into one
tolerance value per position. The question nobody appears to have asked of this dataset is whether
the **identity** of the inserted residue interacts with the site.

**Hypothesis H5.** Which residue is best to insert depends on the local environment: at buried or
structured positions only small and flexible residues are tolerated, while exposed positions accept
anything.

**Test.** For each inserted residue type, the Spearman correlation of its per-position scores with
every other type's; then a two-way analysis of position, inserted residue, and their interaction,
and specifically the correlation between a residue's advantage (its score minus the position's
median) and the position's relative solvent accessibility and pLDDT.

**Marks.** H5 is supported if the interaction is both statistically clear and large enough to
matter: the median pairwise Spearman between inserted-residue types must be below 0.8, **and** the
effect of residue volume on the per-position advantage must differ in sign or by at least 0.2 in
Spearman between buried and exposed positions, split at the median accessibility.

**Why it is worth asking.** If the answer is yes, it is a design rule for what a tag's flanking
linker should be made of, which is the practical question behind placing GFP11, and it comes from
data that is already public and already paid for.

**Limits fixed in advance.** One protein and one assay, so any rule found here is a hypothesis for
other proteins, not a general law. Codon-level effects (translation, mRNA structure) are not
separable from amino-acid effects in this dataset.

## Gate 3 result (25 Nov 2026): partially supported, and Gate 5 refuted

`analysis/07_gate3.py`, output `results/gate3.json`.

**Gate 3, the pLDDT claim on proteins that did not generate it.**

| Protein (assay) | pLDDT | ProDomino | best feature there |
|---|---|---|---|
| Kir2.1, 1-residue (surface expression) - generated the claim | 0.850 | 0.747 | pLDDT |
| **TEM-1** (enzyme fitness), 274 positions | **0.871** | 0.768 | pLDDT |
| TEM-1 excluding the signal peptide, 253 positions | 0.808 | 0.676 | pLDDT |
| **DHFR** (folding sensor), 185 positions | 0.709 | 0.712 | distance to terminus, 0.825 |

**The pre-registered mark required pLDDT to be at least ProDomino on both held-out proteins. It
wins clearly on TEM-1 and is level on DHFR (0.709 against 0.712), so by the rule written in advance
the claim is supported on one of two and is recorded as partial, not as support.**

The weaker statement that does hold across all three proteins and all insert sizes tested: pLDDT is
never worse than ProDomino by more than 0.003 AUROC, and is ahead by 0.10 on two of three. A free
number printed by AlphaFold matches or beats a 3-billion-parameter model on every dataset here.

**The pattern that is now the most robust result of this project.** Across seven proteins, which
feature predicts best keeps changing: conservation on AraC (0.828), distance to the terminus on Flp
(0.882) and on DHFR (0.825 to 0.842 across its four conditions), accessibility on TVMV (0.745),
pLDDT on Kir2.1 (0.850) and TEM-1 (0.871). On DHFR, a 187-residue single domain, simply being near
either end beats every structural and learned predictor. No feature, and no model, is best
everywhere, and a feature that is strong in one protein can be at chance or reversed in another.

**Gate 5, the identity of the inserted residue: refuted.** Over the 4,508 usable TEM-1 insertions,
the median Spearman correlation between the per-position profiles of different inserted residues is
**0.917** (range 0.779 to 0.986). The mark required below 0.8. **Where you insert dominates what you
insert**, by a wide margin, so the interaction hypothesis fails on its first condition and the
second was not evaluated.

A small real effect does sit underneath that, in the expected direction: averaged over all
positions, the most tolerated single insertions are Thr 0.188, Ser 0.178, Ala 0.174 and Gln 0.167,
and the least tolerated are Trp 0.112, Pro 0.119 and Leu 0.126. Small polar residues are best and
tryptophan and proline are worst, which is the textbook expectation and is a weak effect next to
position. For designing a tag linker this says the composition matters much less than the site.

## Gate 4 result (25 Nov 2026): the segment hypothesis is refuted, including the part that first looked right

`analysis/08_gate4_segments.py` and `analysis/p3_check.sh`, output `results/gate4.json`.

**P1 fails, on every dataset.** Run length is a *worse* predictor than the residue's own pLDDT
everywhere, not better: the advantage of `run_len` over pLDDT is -0.065, -0.063 and -0.067 AUROC
for Kir2.1 at 1, 2 and 3 residues, -0.068 on TEM-1, -0.083 on DHFR, and -0.117, -0.062 on AraC and
TVMV, with only Flp (+0.014) and SigF (+0.039) marginally positive. Collapsing a position into the
run that contains it destroys information rather than adding any. The same holds at pLDDT
thresholds set at the 40th and 60th percentile instead of the median.

**P2, the real test, fails.** The advantage does not grow with insert size: within Kir2.1 it is
flat at about -0.065 for 1, 2 and 3 residues. The domain-insertion side nominally passed its clause,
but only because two of its four proteins have small positive values and the comparison is between
different proteins in different assays, so it carries no weight on its own.

**P3 passed and then did not survive a check I ran before reporting it.** Raw Spearman between
slack and size sensitivity is -0.302 (95% -0.382 to -0.220), which meets the mark. But slack and
pLDDT correlate at -0.855, and the partial Spearman of slack against size sensitivity **controlling
for pLDDT is -0.068, 95% interval -0.167 to +0.032, which includes zero**. The apparent effect was
pLDDT in another coordinate system.

**Verdict: H4 is refuted.** Insertion tolerance in these data is described by a property of the
single position, and the length and spare room of the flexible segment around it add nothing once
that position's own pLDDT is known. The physical intuition that an insert needs room to go is not
visible here, at least for inserts of one to three residues; whether it appears for a 16-residue
tag cannot be answered, because no dataset measures different insert sizes at matched positions
above three residues.

## Where the project stands after five gates

| Gate | Question | Result |
|---|---|---|
| 1 | Does ProDomino beat transferable baselines on unseen proteins? | Failed its mark; 0.65-0.67 on unseen against 0.877 on its own protein, and uncalibrated across proteins |
| 2 | Does tolerance depend on insert size? Does a domain model transfer to short inserts? | Transfers (0.747); size matters in the tail (worst-quarter overlap 0.507); pLDDT beat it |
| 3 | Does the pLDDT result hold on proteins that did not generate it? | Partial: wins on TEM-1, ties on DHFR, never worse than -0.003 |
| 4 | Is tolerance a segment property? | **Refuted**, including after controlling the one part that passed |
| 5 | Does the identity of the inserted residue matter? | **Refuted**: median correlation between inserted residues 0.917; where beats what |

**What is solid enough to write down.** Seven proteins, three assay types, two published models'
worth of comparison: no single feature or model predicts insertion tolerance well across proteins,
the best feature changes from protein to protein, a free AlphaFold confidence score matches or beats
the published sequence model everywhere tested, and two mechanistic explanations for the variation
(segments, insert identity) are ruled out.

**What is not there.** A mechanism. The project has found that the field's models are weaker than
advertised and has eliminated two candidate explanations, which is useful and publishable as a
benchmark with negative results, but it has not produced a positive discovery.

## Gate 6: how many measured positions does one protein need? (written 25 Nov 2026, before computing)

Five gates established that no feature and no published model transfers well between proteins. That
is a negative result about general prediction. It implies a different, practical question, which is
the one a laboratory actually faces: **if you must measure some positions in your own protein, how
few can you get away with, and which ones should they be?**

**Design.** For each dataset separately (Kir2.1 at 1 and at 3 residues, TEM-1, DHFR, and the four
DI_screen domain-insertion proteins), a logistic regression on the per-position features already
computed (pLDDT, relative accessibility, loop, distance to terminus, ProDomino) is fitted on k
positions of that protein and scored by AUROC on all its remaining positions. k runs over 8, 16,
24, 32, 48, 64, 96 and 128, with 50 random repeats at each k.

**Selection strategies, all chosen before any measurement**, which is the real constraint - a
laboratory picks its panel first and measures once:

- **random**;
- **stratified**: equal numbers from each decile of pLDDT, so the panel spans the range;
- **diverse**: greedy maximum-minimum distance in standardised feature space.

Sequential active learning is deliberately excluded, because it requires several rounds of
measurement and is not what a laboratory will do.

**Marks fixed now.**

- **M1.** At k = 24, the fitted per-protein model beats the best zero-measurement predictor for that
  protein (pLDDT alone) on at least **5 of the 8 datasets**, judged on medians over the 50 repeats.
- **M2, the interesting one.** At k = 24, the better of stratified and diverse beats random by at
  least **0.02 AUROC averaged over datasets**, with a bootstrap interval over repeats excluding
  zero. If this fails, the honest recommendation is "measure any 24 positions", which is a weaker
  but still usable result.

**Reported without a mark.** The full curve of AUROC against k per dataset; the k at which each
dataset first reaches AUROC 0.85 and 0.90, if it does; which features carry the per-protein models;
and whether including ProDomino as a feature helps once real measurements are available.

**Limits fixed in advance.** Positions within one protein are not independent - neighbours share
structure - so a random split flatters any method that smooths locally; a blocked split by
contiguous stretches is therefore also reported, and is the more honest number. Eight datasets
from five laboratories and four assay types is a small basis for a protocol recommendation.

## Gate 7 and 8, both cheap, run alongside

**Gate 7, consensus.** Since the best single feature changes between proteins, the mean rank of
pLDDT, accessibility, distance to terminus and ProDomino is tested as a zero-training predictor.
Mark: its worst AUROC across the eight datasets is higher than the worst AUROC of each single
feature, i.e. it is more robust even if not best anywhere.

**Gate 8, protein language model.** ESM-2 3B is already loaded for ProDomino. Two further
zero-training predictors are scored: the entropy of the model's residue distribution at each
position, and the log-likelihood of the wild-type residue there, both from one forward pass on the
unmasked sequence. Mark: either beats pLDDT on at least 5 of the 8 datasets. This asks whether a
general protein language model already contains what ProDomino was trained to extract.

## Gate 6 and 7 results (25 Nov 2026): both marks fail, and the failure points somewhere specific

`analysis/10_gate6_how_many.py`, output `results/gate6.json`. Eight datasets, 50 repeats per
point, logistic regression on pLDDT, accessibility, loop, distance to terminus and ProDomino.

**M1 fails: 4 of 8, the mark was 5.** Fitting a model on 24 measured positions of a protein beats
using pLDDT alone on DHFR (0.867 against 0.709), SigF (0.614 against 0.460), Flp (0.830 against
0.749) and AraC (0.828 against 0.779), and loses on Kir2.1 at both insert sizes, TEM-1 and TVMV.

**The pattern is not random: measuring helps exactly where the free predictor is bad.** Every
dataset where the fitted model wins is one where pLDDT alone is weak (0.46 to 0.78); every dataset
where it loses is one where pLDDT is already strong (0.68 to 0.89). On Kir2.1 the fitted model does
not catch pLDDT even with 128 measured positions (0.812 against 0.850).

**The likely reason, and it suggests a better protocol.** Fitting five coefficients on 24 points is
noisy, and the noise costs more than the per-protein adaptation gains when one feature is already
excellent. The protein-specific thing worth learning from a small panel is probably not a set of
weights but **which predictor to trust** - which is the question this whole project has been
circling, since the best feature changes from protein to protein. Choosing among five candidates
needs far less data than fitting five weights. That is Gate 9 below.

**M2 fails: mean gain 0.0116, the mark was 0.02.** Designing the panel (stratified by pLDDT, or
maximum-coverage in feature space) beats picking 24 positions at random by about one AUROC point on
average, positive in six of eight datasets but small: per-dataset gains 0.001, 0.009, 0.004, -0.001,
0.018, 0.056, -0.009, 0.015. The honest recommendation from this is that panel design is worth
little and **any 24 positions will do about as well as a cleverly chosen 24**.

**Gate 7, consensus: the mark fails, and the numbers are still the most useful table in the
project.** Across the eight datasets, as zero-measurement predictors:

| Predictor | mean AUROC | worst case |
|---|---|---|
| **consensus rank of four features** | **0.764** | 0.492 |
| pLDDT | 0.749 | 0.460 |
| **ProDomino** | 0.725 | **0.646** |
| accessibility | 0.692 | 0.481 |
| distance to terminus | 0.656 | 0.351 |
| in a loop | 0.535 | 0.493 |

The mark required the consensus to have a higher floor than every single predictor; ProDomino's
floor is higher, so it fails.

**And that comparison is worth stating fairly, because it cuts against the earlier result.**
ProDomino has the lowest peak of the serious predictors but **the highest floor by a wide margin**:
it never drops below 0.646, while pLDDT falls to 0.460 and distance-to-terminus to 0.351, both
below chance, on SigF. For someone holding an unfamiliar protein and no measurements, a predictor
that is never terrible may be worth more than one that is sometimes excellent. The earlier finding
that a free number matches or beats the model still holds on average and on peak; it does not hold
on the worst case, and that was not tested until now.

SigF is the dataset where everything fails: pLDDT 0.460, accessibility 0.481, loop 0.493, distance
0.351, consensus 0.492, and only ProDomino is above chance at 0.673.

## Gate 9: use a small panel to CHOOSE a predictor, not to fit one (written 25 Nov 2026, before computing)

Gate 6 showed that fitting weights on 24 measured positions loses to a single good feature whenever
that feature is already strong, and that the datasets where fitting wins are exactly those where
the free predictor is weak. Gates 1 to 7 showed that which predictor is best changes from protein
to protein and that nobody can tell in advance which one it will be. Those two facts together
suggest the protein-specific quantity worth learning from a small experiment is **which predictor
to trust**, not a set of coefficients. Choosing among five candidates costs far less data than
fitting five weights.

**This is a post-hoc idea**, generated by Gate 6's failure pattern, and is marked as such. It is
pre-registered here before being computed, but it did not exist before seeing those results, so it
needs replication on datasets not used to form it before it should be believed.

**Procedure.** For each dataset and each of 50 repeats: draw k positions at random, compute each
candidate predictor's AUROC **on those k positions only**, take the best, and report that
predictor's AUROC **on the remaining positions**. k runs over 8, 16, 24, 32, 48, 64, 96, 128.

**Candidates**, fixed now: pLDDT, relative accessibility, distance to terminus, ProDomino, and the
consensus rank of those four. The loop indicator is excluded because Gate 7 showed it is at chance
everywhere; including a known-useless candidate would only add selection noise. ESM-2 entropy and
wild-type log-likelihood are added if Gate 8 shows either is competitive, and that decision is
recorded before the selection results are examined.

**Compared against**, all on the same held-out positions: always pLDDT; always ProDomino; always
consensus; the fitted logistic model of Gate 6; and an oracle that picks the best predictor using
the held-out positions themselves, as an upper bound.

**Marks at k = 24.**

- **M3.** Selection's mean AUROC across the eight datasets exceeds always-pLDDT's (0.749), **and**
  its worst case across datasets exceeds pLDDT's worst case (0.460).
- **M4, the strong one.** Selection's worst case also reaches or exceeds ProDomino's worst case
  (0.646), which would make it the first predictor in this project that is both good on average and
  never bad.

**Reported without a mark.** The k at which selection overtakes always-pLDDT; how often it picks
each candidate per dataset; and the gap to the oracle, which measures how much of the loss is
selection noise rather than a missing predictor.

**Limits fixed in advance.** Selecting on the same small sample that was drawn from the dataset
being evaluated is optimistic relative to a real laboratory, which would measure a panel once with
assay noise this analysis does not model. Eight datasets is a small basis, and three of them
(Kir2.1 at two sizes, and its deletions) are not independent of each other.

## Gate 8 result (25 Nov 2026): ESM-2 features

Computed in `analysis/09_esm_features.py` from one forward pass of ESM-2 3B per protein, as the
entropy of the model's residue distribution and the log-likelihood of the wild-type residue.
Results are reported with Gate 9 below.

## Gate 8 and 9 results (25 Nov 2026): the first positive result of the project

`analysis/11_gate9_select.py` and `analysis/blocked_check.sh`, output `results/gate9.json`.

**Gate 8 fails: a general protein language model is not a good insertion predictor.** ESM-2 3B
entropy averages 0.678 AUROC across the eight datasets and wild-type log-likelihood 0.668, against
pLDDT's 0.748. Entropy beats pLDDT on 2 of 8 datasets and log-likelihood on 1; the mark was 5. It
is better than pLDDT exactly where pLDDT collapses (SigF 0.614 against 0.460, AraC 0.817 against
0.779) and much worse elsewhere. By the rule fixed beforehand, neither was added as a Gate 9
candidate.

**Gate 9 passes both marks, under the random panel and under the harder blocked one.**
At k = 24 measured positions:

| Method | random panel: mean / worst | blocked panel: mean / worst |
|---|---|---|
| **select the predictor using the panel** | **0.816 / 0.666** | **0.796 / 0.654** |
| oracle (upper bound, picks with full knowledge) | 0.825 / 0.670 | 0.828 / 0.672 |
| fit a 5-feature model on the panel (Gate 6) | 0.772 / 0.603 | - |
| always consensus | 0.763 / 0.489 | - |
| always pLDDT | 0.747 / 0.459 | 0.749 / 0.464 |
| always ProDomino | 0.724 / 0.647 | 0.725 / 0.641 |

**M3 passes**: selection beats always-pLDDT on the mean (0.816 against 0.747) and on the worst case
(0.666 against 0.459). **M4 passes**: its worst case also reaches ProDomino's floor (0.666 against
0.647). It is the first predictor here that is both good on average and never bad, and the first
thing in this project to beat every fixed alternative on both axes at once.

**Selection costs almost nothing against the oracle: 0.816 against 0.825.** Twenty-four positions
are enough to identify which predictor to trust; almost all of the remaining gap is a missing
predictor, not selection noise. It overtakes always-pLDDT already at **k = 8**.

**And it recovers the right answer per protein**, which is the clearest confirmation that the
per-protein variation found in Gates 1 to 3 is real signal and not noise. Out of 50 repeats at
k = 24, the predictor chosen most often is:

| Dataset | chosen | times | matches the best feature found earlier |
|---|---|---|---|
| Kir2.1, 1 and 3 residues | pLDDT | 34, 36 | yes |
| TEM-1 | pLDDT | 28 | yes |
| DHFR | distance to terminus | 24 | yes |
| AraC | ProDomino | 33 | yes |
| Flp | distance to terminus | 27 | yes |
| SigF | ProDomino | 40 | yes |
| TVMV | accessibility | 25 | yes |

**Caveats, which matter here more than anywhere else in this file.**

1. **This hypothesis is post-hoc.** It was generated by Gate 6's failure pattern and tested on the
   same eight datasets. It needs replication on datasets not used to form it before it should be
   believed, and that is the single most important next step.
2. The blocked check weakens but does not remove the result (mean 0.816 to 0.796), and its margin
   over ProDomino's floor becomes thin (0.654 against 0.641). Flp contributes only 9 usable repeats
   under blocking, because contiguous stretches of it often contain none of its 30 positives.
3. No assay noise is modelled: a real panel of 24 measurements carries error this analysis does not.

**What it means in one sentence.** There is no universal rule for where a protein tolerates an
insertion, and this project spent five gates confirming that - but you do not need one: about two
dozen measured positions are enough to find out which of the existing rules applies to your
protein, and that beats trusting any of them in advance.

# Replication, 26 Nov 2026

## Gate 10: replicate Gate 9 on the Tsuboyama set (written before any score was computed)

Gate 9 was post-hoc and was tested on the same eight datasets that produced it. This replicates it
on data that played no part in forming it.

**The replication set.** Tsuboyama et al. 2023 (Zenodo 7992926), library 3. Each row is one variant
of one small domain with a measured folding free energy. Extracted with
`analysis/13_tsuboyama_extract.py` and prepared with `analysis/14_tsub_prepare.py`:

- **401 domains**, 22,099 single glycine insertions and 22,037 single alanine insertions, plus
  18,000 single-residue deletions;
- folding free energies with a median 95% interval width of 0.144 kcal/mol;
- the authors' own AlphaFold model for every domain, and **the numbering check passes perfectly**:
  every domain's deletion rows name the wild-type residue at each position, and all 401 match their
  model at offset zero with no mismatches;
- wild-type dG recovered for 387 of the 401 domains from all four libraries.

This is roughly fifty times the data Gate 9 was tested on, from a different laboratory, a different
assay (proteolysis-based folding stability rather than surface expression, enzyme fitness or
transcription-factor function), and different proteins.

**Label, fixed now.** Tolerance is ddG = dG(variant) - dG(wild type) for the same domain, and a
position is **tolerated if ddG > -1 kcal/mol**, the conventional threshold for a destabilising
mutation. Only domains with a wild-type dG are used. Variants whose dG 95% interval is wider than
1.0 kcal/mol are dropped as unreliable. A domain enters the analysis only if it keeps at least 20
scored positions and has both classes present.

**Primary test: Gate 9's procedure, unchanged.** For each domain and each of 50 repeats, draw k
positions, compute each candidate's AUROC on those k only, take the best, and score that candidate
on the domain's remaining positions. Candidates are exactly those of Gate 9: pLDDT, relative
solvent accessibility, distance to the nearest terminus, ProDomino, and their consensus rank.
k = 24 is the headline, with the full curve reported.

**Marks, the same as Gate 9's and fixed before computing.**

- **R1.** Selection's mean AUROC across domains exceeds always-pLDDT's, **and** selection's mean
  exceeds always-ProDomino's.
- **R2.** Selection's lower decile across domains exceeds always-pLDDT's lower decile. The worst
  case of a single domain is too noisy with only about 37 positions each, so the tenth percentile
  across domains replaces Gate 9's outright worst case; this substitution is recorded before
  seeing any result.

**Also tested, each with its own mark.**

- **R3, the per-protein variation itself.** Gates 1 to 3 found that the best predictor changes from
  protein to protein. Mark: across these domains no single candidate is the best for more than 60%
  of them. A single candidate winning nearly everywhere would mean the earlier variation was an
  artefact of comparing four laboratories rather than a real property.
- **R4, Gate 5 again, on 397 domains instead of one.** Gate 5 found that where you insert matters
  far more than what you insert, from a median correlation of 0.917 between inserted residues in
  TEM-1. Here glycine and alanine insertions are measured at the same positions in 397 domains.
  Mark: the median across domains of the per-domain Spearman between the glycine and alanine
  profiles is above 0.8, which would replicate Gate 5's refutation on a far larger sample.

**Limits fixed in advance.** These domains are small, median length well under 100 residues, and
selected for being stable enough to measure, so they are not a random sample of proteins. All
inserts are one residue. Folding stability is not the same endpoint as function or expression, so
agreement with the earlier gates is evidence that the selection procedure transfers across
endpoints, while disagreement might be about endpoints rather than about the procedure.

## Gate 10 result (26 Nov 2026): the replication fails, and explains why Gate 9 looked good

`analysis/15_gate10_replicate.py`, output `results/gate10.json`. 355 domains, 19,117 positions,
median 57 scored positions per domain.

| Mark | Result | |
|---|---|---|
| R1, selection beats always-pLDDT and always-ProDomino on the mean | 0.859 against 0.852 and 0.783 | passed, by 0.007 |
| **R2, selection's lower decile beats always-pLDDT's** | **0.698 against 0.709** | **failed** |
| R3, no candidate is best in more than 60% of domains | pLDDT best in 164 of 355, 46.2% | passed |
| R4, median Gly-against-Ala agreement above 0.8 | 0.856 over 374 domains | passed |

**Both R1 and R2 were required, so the replication fails.** The gain over simply always using
pLDDT is 0.007 AUROC and the lower decile is worse, which is not a reason to run 24 experiments.

**The k curve shows what happened, and it is the useful part:**

| k | selection | always pLDDT | oracle | gain | headroom |
|---|---|---|---|---|---|
| 8 | 0.819 | 0.845 | 0.894 | **-0.026** | 0.049 |
| 12 | 0.839 | 0.846 | 0.897 | -0.007 | 0.051 |
| 16 | 0.848 | 0.846 | 0.899 | +0.002 | 0.053 |
| 24 | 0.859 | 0.852 | 0.908 | +0.007 | 0.056 |
| 32 | 0.862 | 0.855 | 0.915 | +0.007 | 0.060 |

An oracle that always picks the right predictor gains 0.056 AUROC over always-pLDDT, so the
per-protein variation is real and worth about that much. **A 24-position panel recovers 0.007 of
it, about one eighth.** Below 16 positions, selecting is worse than not selecting at all.

**Why Gate 9 looked better.** Its eight datasets had 185 to 434 scored positions each; these
domains have a median of 57. A 24-position panel is a twentieth of a large dataset and nearly half
of a small domain, so both the panel and what remains to evaluate on are noisy here. Gate 9's
result was not wrong about its datasets, but it does not generalise to small domains, and the
honest recommendation for a small domain is to skip the panel and use pLDDT.

**What does replicate, and is now the firmest result in this project.** On 355 domains from an
independent laboratory, with folding stability rather than expression or activity as the endpoint:

| Predictor | mean AUROC | lower decile |
|---|---|---|
| **pLDDT** | **0.845** | 0.708 |
| consensus of four | 0.834 | 0.670 |
| **ProDomino** | 0.775 | 0.565 |
| distance to terminus | 0.720 | 0.442 |
| relative accessibility | 0.672 | 0.551 |

A number AlphaFold already prints beats the purpose-built 3-billion-parameter insertion model by
0.070 AUROC on the mean and 0.143 on the lower decile. Combined with Gates 2, 3 and 7, this now
holds across eleven proteins or domain sets, four laboratories and four endpoints.

**R4 replicates Gate 5 at scale.** Glycine and alanine inserted at the same positions in 374
domains agree with median Spearman 0.856 (quartiles 0.770 and 0.906), and mean ddG differs by only
0.12 kcal/mol (-2.03 for Gly, -1.91 for Ala). Where you insert dominates what you insert, now on
374 domains instead of one protein.

## Gate 11: does the preferred inserted residue depend on local structure? (written 26 Nov 2026, before computing)

**Literature checked first this time**, at the user's instruction, and this question builds on it
rather than rediscovering it.

What is already known:
- Gutierrez & Rocklin (bioRxiv 2024.12.18.629072) analysed the **103 stabilizing indels** of this
  same dataset and did compare glycine with alanine, but only at the **36 sites that an insertion
  stabilises**. They report "Gly is nearly always tolerated at sites that are stabilized by an Ala
  insertion" while "Ala is often destabilizing at sites stabilized by a Gly insertion", and explain
  it by glycine's wider backbone range, with an SH3 case where glycine relieves strain next to a
  positive-phi residue. They state their analysis is confined to the stabilizing subset.
- Helix propensity is textbook: alanine highest, glycine lowest among non-proline residues, about
  1 kcal/mol apart, with glycine favoured at helix caps (Serrano & Fersht 1992; Pace & Scholtz 1998).
- Topolska et al. 2025 report that insertions and deletions inside helices are roughly three times
  as destabilising as substitutions there.

**What is not done**: a position-matched glycine-against-alanine comparison across all positions,
resolved by local secondary structure. This dataset has both residues inserted at the same ~22,000
positions in 374 domains, and Gate 10 already shows the two agree only moderately (median Spearman
0.856), so there is per-position disagreement to explain.

**A tension worth testing.** Gutierrez & Rocklin find glycine more permissive at stabilizing sites,
yet across all positions here alanine is the milder insertion on average (mean ddG -1.91 against
-2.03 kcal/mol). Both can hold only if the preferred residue **changes with context**.

**Hypothesis H11.** The preferred inserted residue depends on local secondary structure: alanine is
relatively favoured inside helices, glycine in loops and turns.

**Measure.** For every position with both measurements, delta = ddG(Ala) - ddG(Gly); positive means
alanine is the milder insertion there. Secondary structure of the position is taken from the
authors' AlphaFold model by backbone dihedrals, three classes: helix, strand, loop.

**Marks, fixed now.**
- **A.** Mean delta in helix positions exceeds mean delta in loop positions, with a 95% interval
  from bootstrapping **over domains** that excludes zero. Bootstrapping over domains, not
  positions, because positions within a domain are not independent.
- **B.** The difference between the helix and loop means is at least **0.3 kcal/mol**, which is a
  meaningful fraction of the roughly 1 kcal/mol helix propensity gap and large enough to matter
  when choosing a tag linker.

Both are required. If A holds and B does not, the effect is real but too small to act on, and that
is what will be reported.

**Secondary, reported without a mark.** Delta against burial; delta at positions whose flanking
residues sit in the positive-phi region, where glycine should be favoured most sharply if the
Gutierrez & Rocklin mechanism generalises; and whether the helix effect concentrates in the middle
of helices rather than at their ends.

**Limits fixed in advance.** Small stable domains, folding stability as the endpoint, single-residue
inserts of only two chemical types, and secondary structure assigned from a predicted model by
dihedrals rather than DSSP on an experimental structure. A positive result is a statement about
these domains and this endpoint, and the practical reading for tag design would still need testing
with a real tag.

## Gate 11 result (26 Nov 2026): the direction is confirmed, the pre-registered effect size is not

`analysis/16_gate11_gly_ala.py` and `analysis/g11_control.sh`, output `results/gate11.json`.
19,291 position-matched glycine/alanine pairs in 385 domains. delta = ddG(Ala) - ddG(Gly), so
positive means alanine is the milder insertion at that position.

| Local structure | n | mean delta (kcal/mol) |
|---|---|---|
| helix | 9,312 | **+0.201** |
| strand | 7,930 | +0.102 |
| **loop** | 2,049 | **-0.082** |

**Mark A passes.** Helix minus loop is +0.283 kcal/mol, 95% interval by bootstrapping over domains
[0.244, 0.321], which excludes zero. The sign flips between helix and loop, as predicted.

**Mark B fails.** The pre-registered threshold was 0.3 kcal/mol and the result is 0.283. Both marks
were required, so **H11 as written is not supported**, by 0.017 kcal/mol.

**The three component effects are real and independent of each other**, which was the check that
could have killed all of this:

- **Helix against loop holds inside every burial tertile**: buried +0.356, exposed +0.221, middle
  +0.199. It is not a burial artefact.
- **Burial holds inside every structure class**: strand +0.129, helix +0.197, loop +0.062. Alanine
  is favoured when the position is buried, which is what a methyl group filling hydrophobic space
  should do.
- **The positive-phi effect holds inside every structure class**: where a flanking residue sits in
  the left-handed region, glycine gains by 0.310 in strands, 0.259 in helices and 0.158 in loops
  (overall difference -0.294, interval [-0.327, -0.262], 2,163 such positions).

**The clearest post-hoc contrast**, reported as post-hoc: helix middles against loops is **+0.381
kcal/mol**, interval [0.339, 0.423]. Alanine's advantage sits almost entirely in the middle of
helices (+0.299) and vanishes at helix ends (+0.000), which is what the classical helix propensity
work predicts, since glycine is tolerated at helix caps.

**What this adds to what was already known.**

1. Helix propensity, established for substitutions, carries over to the identity of an **inserted**
   residue, and it is localised to helix middles.
2. Gutierrez & Rocklin proposed from one SH3 case that glycine wins where the backbone needs a
   left-handed conformation. That generalises: it holds over 2,163 positions and, unexpectedly,
   **inside strands and helices as strongly as in loops**, so it is not a loop phenomenon.
3. Which residue to insert is a secondary consideration. Gate 10 showed position dominates identity
   (median agreement 0.856), and the whole range of this effect, from +0.30 in buried helix middles
   to -0.22 in loops with a left-handed neighbour, is about 0.5 kcal/mol against a mean insertion
   cost near -2.0 kcal/mol. **Choosing the site is worth about four times as much as choosing the
   residue.**

**Practical reading, with that caveat.** When a position must be used and the residue can be
chosen, use alanine inside a helix or at a buried position, and glycine in a loop or next to a
residue in the left-handed region. The gain is real and modest.
