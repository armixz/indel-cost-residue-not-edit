# Is it worth building a tag-insertion data company? Research and one decisive test

Date: 2026-11-26. Pre-registration: `PREREG_abundance_proxy.md` (written before any data
was retrieved). Analysis: `analysis/46_abundance_proxy.py`, `analysis/47_indel_proxy.py`.
Raw output: `results/abundance_proxy.json`, `results/indel_proxy.json`.

**Answer: no. Not as a data company. The reason is not cost and not competition; it is
that the thing that would make the data cheap to produce does not work.**

---

## 1. The question this reduced to

Producing library-scale data on where a peptide tag can be inserted inside a protein is
cheap except for one step: measuring whether the host protein still works. A bespoke
activity assay costs 1-3 months per protein, and that is the whole reason the dataset does
not exist.

Every general, protein-agnostic, pooled assay that exists measures **abundance or
stability**, never activity: VAMP-seq, cDNA-display proteolysis, degron reporters,
split-GFP, HiBiT. So the business case rests on a single empirical claim:

> if a variant is still abundant, it is still functional.

If that holds, one general assay replaces the bespoke ones and the idea scales. If it does
not, there is no shortcut, and the expensive step is exactly the step that cannot be
amortised across customers.

## 2. The test

Four proteins have an abundance deep mutational scan and an independent activity scan over
the same variants, all from MaveDB. Scores were put on a common 0 = null, 1 = wild-type
scale; nonsense variants were dropped before scoring, because they are trivially both
non-abundant and non-functional and would inflate every correlation.

| protein | n (missense) | Spearman rho | AUC | false-pass rate | % non-functional |
|---|---|---|---|---|---|
| PTEN | 4,839 | 0.325 | 0.709 | 14.2% | 24.9% |
| NUDT15 | 2,844 | 0.340 | 0.699 | 14.3% | 23.2% |
| CYP2C9 | 4,421 | 0.749 | 0.891 | 16.3% | 46.4% |
| G6PD | 8,217 | 0.373 | 0.704 | 30.5% | 46.1% |

**Median rho 0.357. Only 1 of 4 proteins reaches the pre-registered pass bar of 0.7.**

The **false-pass rate** is the number that matters: among variants an abundance assay would
call fine (normalised abundance >= 0.75), the fraction that are actually non-functional
(normalised function < 0.5). Pooled, **22.3%**.

Pre-registered interpretation, applied mechanically: **FAIL**.

### Threshold sensitivity

The cutoffs were chosen in advance, so the result is only worth something if it survives
other reasonable ones. False-pass rate over the grid (pooled, 20,321 variants):

| abundance cut | func<0.3 | func<0.4 | func<0.5 | func<0.6 | func<0.7 |
|---|---|---|---|---|---|
| >= 0.60 | 0.169 | 0.202 | 0.241 | 0.292 | 0.362 |
| >= 0.70 | 0.161 | 0.191 | 0.228 | 0.277 | 0.345 |
| >= 0.75 | 0.159 | 0.187 | 0.223 | 0.271 | 0.337 |
| >= 0.80 | 0.157 | 0.185 | 0.219 | 0.266 | 0.330 |
| >= 0.90 | 0.166 | 0.191 | 0.221 | 0.263 | 0.320 |

The best cell on the entire grid is **15.7%**. There is no choice of thresholds that
rescues abundance as a proxy.

The more damaging detail is the near-flat columns: **raising the abundance bar from 0.60
to 0.90 changes the false-pass rate from 0.169 to 0.166.** Demanding that a variant be
*more* abundant buys essentially no additional confidence that it still works. That is
what it looks like when two measurements are weakly related rather than noisily related,
and it rules out the obvious fix of running the cheap assay at higher stringency.

### Direction of the bias

A substitution is a smaller perturbation than an 11-16 residue tag. This asymmetry was
written into the pre-registration before any result was seen: a bad false-pass rate on
substitutions is bad news that carries over, while a good one would only have been weak
good news. The result came out bad, so it carries.

## 3. Two datasets retrieved and then rejected, with reasons

Both rejections are recorded because each one first produced a result that looked real.

**ASPA** (abundance `00000657-a-1`, "function" `00000657-b-1`). Its second assay is a
cellular **toxicity** readout, not an activity readout. In a toxicity assay a premature
stop makes no protein and is therefore not toxic, so nonsense variants score as benign
(median 0.122) while missense score lower (0.048). A null allele lands on the good end,
the assay has no loss-of-function orientation, and it cannot answer this question.
Direction of the cut, recorded so it cannot be read as convenient: including ASPA with its
sign flipped would have contributed rho = +0.77, the second highest of the set, and raised
the median. **Dropping it works against the conclusion reached here, not for it.**

**GRB2** (`00001305-a-1` / `-b-1`, LABEL-seq). This one was pursued hard, because its
library contains 191 single-residue deletions and 119 delins variants with *both* assays
on the same library -- a direct test on length-changing variants rather than substitutions.
It failed three ways: (a) **0.0%** of variants fall below the loss-of-function threshold,
so the activity assay has essentially no dynamic range for loss of function; (b) nonsense
variants score at or above wild type in *both* assays (abundance 1.014 vs synonymous
0.990), meaning the abundance readout does not detect truncation at all, consistent with
the signal coming from a terminal tag a premature stop leaves intact; (c) abundance and
activity are **anti**-correlated throughout (rho about -0.2 including substitutions), the
signature of an overexpression system where more GRB2 sequesters signalling components and
suppresses ERK.

A normalisation bug nearly turned (b) into a finding. Anchoring on nonsense variants, as
script 46 does, gave a denominator of (0.990 - 1.014) = **-0.024** -- tiny and negative --
which inverted and exploded the scale and made every correlation in the first run
negative. That was an artefact of my own rescaling, caught by checking the anchors against
the raw class medians rather than by anything in the output looking wrong.

**Consequence: the conclusion rests on four proteins and substitutions. There is still no
public dataset that measures abundance and function for multi-residue insertions in the
same library.** That gap is real, and it is the honest limit of this test.

## 3b. Independent confirmation, and a larger published version of the same test

A literature search run in parallel with this analysis, and blind to it, turned up the
published value for one of the four pairs: Amorosi et al. (Am J Hum Genet 2021) report
CYP2C9 abundance vs activity at **Spearman 0.749 over n = 4,421**. This analysis computed
**0.749 over n = 4,421** from the raw MaveDB score sets, without having seen that number.
An exact match on an independently derived figure is the best evidence available that the
pipeline here is doing what it claims.

The same search found that a much larger version of this test already exists and is more
pessimistic than mine: Gerasimavicius, Livesey & Marsh (Protein Science 32:e4688, 2023)
compared predicted stability effects against functional DMS scores over **45 datasets, 35
proteins, 161,441 variants** and report a mean Spearman of **0.26** — below my 0.357, and
*worse* than a generic variant-effect predictor at around 0.35.

Three further published figures agree on the asymmetry this analysis found, each from a
different protein and assay:

| finding | number | source |
|---|---|---|
| TYK2 signalling-dead variants with no abundance decrease | **948/2,578 = 37%** | Howard et al., bioRxiv 2025 |
| PTEN/NUDT15 broken variants invisible to abundance | **40-50%** | Cagiada et al., Mol Biol Evol 2021 |
| PTEN clinical discrimination, abundance alone vs multidimensional | **AUROC 0.62 vs 0.92** | Pendyala et al., Cell 2026 |
| pathogenic missense explained by destabilisation, 503 domains | 60%, so **40% unexplained** | Beltran et al., Nature 2025 |

A fifth estimate, and methodologically the cleanest of them, comes from a true pooled DMS
carrying both readouts on the same library: Klesmith et al. (PNAS 114:2265, 2017) measured
a solubility proxy (yeast surface display, and separately Tat selection) alongside enzyme
fitness for two enzymes, and the **solubility proxy explained 37% of the variance in the
fitness landscape**. Together with CYP2C9 (R^2 = 0.56), PTEN/NUDT15 (50-60%), TYK2 (63%)
and the Domainome (60%), five independent estimates converge on the same place: **a
folding or abundance proxy captures roughly half of function.**

**The mechanism, demonstrated directly.** Foit et al. (Mol Cell 36:861, 2009) built a
stability reporter by inserting a test protein into TEM-1 beta-lactamase at 196/197;
stability versus conferred resistance gave R^2 = 0.89, so the reporter worked exactly as
designed. The stabilising mutations it then selected **map predominantly to surface
residues in the protein's binding site and destroy its function.** The assay does not
merely miss function-disrupting variants; under selection it actively moves away from
function. This matters here because it is the same failure mode a surface-loop tag
insertion produces, shown in a system mechanistically unrelated to PTEN — so the
abundant-but-dead class is not a PTEN peculiarity.

So the conclusion is not resting on my four proteins. The published literature reaches it
at larger scale, and the canonical counterexamples are structural rather than noise: the
PTEN dominant negatives p.Cys124Ser, p.Arg130Gly and p.Pro38Ser sit **at or above
wild-type abundance** while being catalytically dead.

**One number that does not match, reported because it does not.** Matreyek et al. (Genome
Medicine 2021) classify PTEN variants four ways and report **6%** abundant-but-dead, where
this analysis gets 14.2% for PTEN. The definitions differ -- their classification uses its
own cutoffs and variant set -- and my own threshold sweep shows this quantity moves by a
factor of two across reasonable cutoffs. It is a reminder that the precise percentage is
soft; what is robust is that it never approaches zero.

## 3c. A confound specific to internal insertions, which I had not accounted for

Terminal-fusion abundance assays have never had to separate three different causes of low
abundance, because only one of them could occur. An **internal** tag insertion creates all
three: (1) the host unfolded; (2) the inserted tag or its new junction **is itself a
degron**; (3) the insertion **exposed a pre-existing cryptic degron**.

This is not hypothetical. Global protein stability profiling has found on the order of
**10^5** proteome-derived short peptides that destabilise a fused GFP, and in yeast the
quality-control degrons are **predominantly internal and alpha-helical** rather than
terminal. No published pooled assay separates these three causes; doing so would need
proteasome inhibition or an E3-ligase-knockout background as an orthogonal arm, which has
never been built into a pooled stability assay.

So an abundance readout on an internal tag-insertion library is a *worse* instrument than
the one measured in section 2, not an equally good one. The 22% false-pass rate is a
floor.

## 4. What the literature search added

Four parallel searches. The findings that change the decision:

**Insertion tolerance does not transfer — and this is what kills the data product.**
- Between two different inserted domains in the *same* host: 69/229 sites (30%) are
  differentially permissive (Coyote-Maestas et al., Nat Commun 10:290, 2019).
- Between homologous host proteins: "functional insertion sites for each protein did not
  overlap substantially" (Nadler et al., Nat Commun 7:12266, 2016).
- Tolerance "is not correlated with conservation, surface exposure, or other static,
  structural features" (Hoffmann & Schmidt, Mol Ther Methods Clin Dev 31:101143, 2023).

A data business works by building a corpus once and selling it many times. Here the corpus
for protein A tells a customer holding protein B very little. So the product is not a
dataset; it is a per-target assay run — and the costly part of a per-target run is the
protein-specific functional selection, the one cost that cannot be amortised.

**A contradiction in the literature, flagged because it cuts against the above.** Two of
the parallel searches returned opposite findings on whether a tolerance map transfers
between different *inserted* elements:

- Coyote-Maestas et al. (Nat Commun 10:290, 2019): in Kir2.1, **69/229 sites (30%) are
  differentially permissive** between a PDZ domain and Cib81. Tolerance is cargo-specific.
- Mathony et al. (Adv Sci 10:e2303496, 2023): across AraC, Flp, TVMV protease and SigF,
  insertion sites selective for a single insert domain "were hardly observed" — permissive
  sites are largely **insert-identity-agnostic**.

I cannot adjudicate these from here, and the honest position is that **transfer between
tags is unresolved**, with one study for each answer. If Mathony is right, a map made with
one tag is more reusable than I argued, and the per-tag part of the non-transfer argument
weakens.

What this does *not* touch: transfer between **host proteins**, where both Nadler 2016 and
Coyote-Maestas 2019 agree that sites do not overlap, and the test in section 2, which is
about abundance versus function and is independent of either claim. The recommendation
below rests on those two, not on the per-tag question.

**The folding/function gap is directly documented, not inferred.** In Kir2.1, Pro186
accepted every insertion tested and trafficked normally, yet a flexible linker there
abolished channel function. Surface trafficking is not a proxy for function, in the one
system where both were measured.

**No general "host still works" readout exists.** Across roughly ten host proteins ever
scanned with a folded domain or tag in a pooled format (MBP, TMBP, SpyCas9, Kir2.1, AraC,
Flp, TVMV protease, SigF, TEM-1, adenylate kinase), the functional readout was built from
scratch for that one protein every single time. The only reused readouts are anti-FLAG
surface labelling, which works only for membrane proteins and reports trafficking, and
"does the virus still replicate", which works only inside a virus. The one attempt at a
general folding channel (Zutz et al., Nat Commun 12:6093, 2021, ibpA-GFP) is stated by its
own authors to be non-comparable between unrelated proteins.

**The niche is open, and the confirmed negative is clean:** nobody has run a pooled screen
inserting GFP11 or HiBiT at many internal positions. HiBiT additionally cannot be the
readout of a pooled FACS screen at all, because luminescence is not sortable.

**The demand side is weaker than the science side.**
- OpenCell tagged 1,757 genes using only N- and C-termini, attributed its 25% failure rate
  to low protein abundance and poor HDR, and never names insertion-site choice as a
  bottleneck (Cho et al., Science 375:eabi6983, 2022).
- Only ~0.66% of reviewed human proteins (135 of 20,429) have both termini chemically
  blocked by a signal peptide and a GPI anchor. When one terminus fails, the field's answer
  is the other terminus, which is free.
- Biosensor labs do not screen for insertion sites. jGCaMP8 screened 813 variants of which
  only 64 were linker constructs; dLight1 picked its site by sequence alignment in an
  afternoon. They guess the site, then burn thousands of experiments tuning around it.
- Hand-picked structure-guided loops already hit 20-40% on 3-6 constructs.
- Every comparable dataset was released free within one publication cycle: Tsuboyama,
  the Lehner indel set, Schmidt's AAV insertional profile, all five Ginkgo GDPa sets, and
  Dyno's own 201,426-variant AAV2 insertion library, which is now a free ProteinGym
  benchmark. Ginkgo does not sell datasets; it gives them away and sells assay runs.
- ProDomino reached Nature Methods using **zero new wet-lab data**, training on 174,872
  natural intradomain insertions.
- US academic market: ~41 NIH-funded sensor projects totalling $26.5M/yr, median award
  $450k/yr. Realistic willingness to pay is a few thousand dollars, grant-cycle gated.

## 5. Verdict

The scientific gap is real and was confirmed from several directions. The business is not.

The test says the cheap general assay cannot stand in for the expensive specific one: at
best 1 in 6 variants it calls fine are dead, and running it more stringently does not help.
The literature says tolerance does not transfer between hosts or between tags, so even
data produced the expensive way does not become a reusable corpus. The demand search says
the flagship user of tagging never hit this problem, and that everyone who produces
comparable data gives it away.

**Recommended: do not pursue this as a data company.** The two things that survive:

1. **The negative result is useful, but it is no longer novel enough to be its own paper.**
   The finding reproduces across four proteins and 20,321 variants and the error does not
   shrink with stringency — but Gerasimavicius 2023 already made the same point over
   161,441 variants, and the 37-50% figures are published for TYK2, PTEN and NUDT15. What
   remains unpublished is narrower: **the threshold-sweep result that stringency does not
   help** (0.169 -> 0.166 from a 0.60 to a 0.90 abundance cut), which is the part that
   actually kills the engineering shortcut, and the specific application to insertion data
   platforms. That is a paragraph in the indel manuscript's discussion, or a short
   commentary, not a paper of its own. Claiming more would be overclaiming.

2. **The assay question, if ever revisited, is the only thing worth building.** Not the
   data. A pooled readout that reports tag-in and host-works for an arbitrary protein would
   be the actual invention; the data would be its by-product. Nothing in the 2024-2026
   literature proposes one, and this test shows why the obvious shortcut fails.

   **One genuinely encouraging result sits exactly in this gap, and it should be recorded
   even though it does not change the verdict.** Xu et al. (PLoS Biology 22:e3002615,
   2024) inserted GFP11 into the H1-S2 loop of tubulin with 16-residue GS linkers, with
   the site chosen from alignment variability plus an AlphaFold2 linker-length prediction.
   Embryonic lethality after tbb-1 RNAi was **5.3% for the internal GFP11 tag versus 2.7%
   untagged** — against **76.1% and 99.8%** for N- and C-terminal full GFP. Five tubulin
   isotypes, three species. So an internal 15-residue tag can preserve function far better
   than the terminal tagging the field defaults to, which is the premise of this whole
   idea and it is correct.

   What it is not: every construct was built individually by CRISPR knock-in, about four
   sites and four linker lengths were tried, it was never pooled, and it was never used as
   a folding readout. Note also that GFP11 with no linker does not complement at all,
   so an internal-tag library has a linker-design dimension that terminal fusions do not —
   a tolerance map would have to be measured per linker length, multiplying the
   experiment.

   And it does not touch the finding above: tubulin function here was measured by an
   organism-level phenotype (embryonic lethality), which is exactly the bespoke,
   non-portable readout that makes this expensive.

## 6. What a reviewer should check

- The four protein pairs are independent studies; the correlations are not inflated by
  shared batch effects, but PTEN's two scans come from the same lab.
- 20,321 variants across four proteins is a narrow base for a general claim about
  abundance assays. Six to ten pairs would be better; those are the only four that exist
  with both readouts and a clean loss-of-function orientation.
- The substitution-to-insertion extrapolation is argued by direction, not measured. The
  measurement does not exist publicly, which is itself part of the finding.
