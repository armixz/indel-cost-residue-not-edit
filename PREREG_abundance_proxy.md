# Pre-registration: is abundance a valid general proxy for function?

Written 2026-11-26, BEFORE looking at any of the data below.

## Why this test decides the question

We are evaluating whether it is feasible to produce library-scale data on where a short
peptide tag can be inserted inside a protein. The identified bottleneck is the functional
readout: a bespoke activity assay per protein costs 1-3 months, which does not scale.

The only escape is a general, protein-agnostic, pooled readout. Every candidate general
readout we know of (VAMP-seq abundance, cDNA-display proteolysis stability, degron
reporters, split-GFP solubility, HiBiT quantitation) measures **abundance or stability**,
not **activity**. So the whole business case reduces to one empirical question:

> If a variant is still abundant/stable, is it still functional?

If abundance is a good proxy, one general assay suffices and the idea scales.
If abundance is a poor proxy, every protein needs its own assay and it does not.

This is answerable from public data, because a handful of proteins have BOTH an
abundance deep mutational scan and an independent activity deep mutational scan over
the same variants.

## Data

MaveDB (https://mavedb.org) score sets. Target pairs sought: a VAMP-seq-style abundance
scan and an activity/function scan on the same protein. Candidates known to exist:
PTEN, NUDT15, TPMT, VKORC1, CYP2C9, ASPA, MSH2, CALM1. Final list is whatever pairs are
actually retrievable; the list of candidates is fixed here before retrieval.

Inclusion: single amino-acid substitutions present in both score sets of a pair.
Exclusion: synonymous and nonsense variants (nonsense variants are trivially both
low-abundance and non-functional and would inflate every correlation; this is the single
most important exclusion and it is declared here in advance).

## Marks

Scores are oriented so that HIGH = more abundant / more functional, and both are
normalised so that the synonymous median is 1 and the nonsense median is 0 where those
classes exist; otherwise min-max on the 1st/99th percentile.

- **A1 (correlation).** Spearman rho between abundance score and function score, per
  protein, over missense variants only. Bootstrap CI over variants.
  - Pass as "abundance is a usable proxy" if rho >= 0.7 in a majority of proteins.
  - Fail if the median rho across proteins is < 0.5.

- **A2 (the decisive asymmetry).** Classify variants as functional vs non-functional at
  a threshold of 0.5 on the normalised function score. Among variants that are
  **abundant** (abundance score >= 0.75, i.e. the assay would call them "fine"), what
  fraction are **non-functional**? This is the false-pass rate: the fraction of insertion
  sites a general abundance assay would wrongly advertise as safe.
  - Report per protein and pooled.
  - A false-pass rate above 20% means an abundance-only product ships wrong answers one
    time in five, which we judge in advance to be commercially fatal for a data product
    sold on its accuracy.

- **A3 (discrimination).** AUC of abundance score for predicting loss of function,
  computed from the Mann-Whitney U statistic. Reported with a bootstrap CI.

- **A4 (does structure context help?).** Within each protein, split variants by relative
  solvent accessibility / burial if structure is available. Prediction: abundance tracks
  function well at buried positions (where loss of function is via destabilisation) and
  poorly at exposed positions (where loss of function is via disrupting an active site or
  interface while the fold survives). If this holds, it tells us precisely WHERE an
  abundance-only assay is trustworthy, which is itself a usable product boundary.
  - This mark is exploratory and is labelled as such in any writeup.

## Pre-declared interpretation

- A1 passes and A2 false-pass < 10%: abundance is a good proxy. The general assay exists,
  the idea scales, proceed.
- A1 fails or A2 false-pass > 20%: abundance is not a substitute for function. Any claim
  that a general assay solves this is wrong, and the idea is a protein-by-protein service
  business, not a data platform.
- Intermediate: the honest answer is "abundance is a screen, not a verdict", and the
  product must be scoped to the regime where A4 says it works.

## Caveat fixed in advance

Substitution data is a proxy for insertion data. A single amino-acid substitution is a
smaller perturbation than inserting 11-15 residues. The direction of that bias is
knowable: insertions are larger perturbations, so if anything they will break function
through non-abundance routes at least as often as substitutions do. Therefore a bad
false-pass rate on substitutions is bad news that carries over; a good false-pass rate on
substitutions is only weak good news. This asymmetry is declared before seeing results so
it cannot be deployed selectively afterwards.
