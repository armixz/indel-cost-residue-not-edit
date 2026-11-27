# Pre-registration 5: does the structural pattern of DELETION cost depend on which residue is deleted?

Written 26 Nov 2026, before running anything in this file.

## Why

The result in RESULTS_indel.md - that the secondary-structure pattern of indel tolerance changes
with the inserted residue - rests on a single contrast, alanine against glycine, because those are
the only two residues inserted anywhere in the Tsuboyama libraries (50,049 Ala and 50,139 Gly
insertions; no others). One contrast cannot characterise a dependence.

The deletion side is the mirror image and is far better powered: the deleted residue is the
wild-type one, so all 20 types are present (E 5,947, K 4,873, L 4,239, V 3,839, A 3,479, I 3,124,
T 3,105, G 2,950, ...). If "the structural pattern depends on the residue" is a real property of
single-residue indels rather than a quirk of alanine, it should appear here too.

**H-del: the cost of deleting a residue from a helix, relative to deleting one from a loop,
depends on that residue's helix propensity; the same for strands and beta-sheet propensity.**

## The confounder that must be handled

The deleted residue is not assigned at random. Evolution put high-helix-propensity residues into
helices, and buried positions hold large hydrophobics. Two separate things could produce an
apparent interaction:

- **burial**: helix interiors and strand cores are more buried, and deleting a buried residue
  costs more whatever it is
- **volume**: strands are rich in beta-branched and large residues

So the test is the **interaction** with structure after those main effects, not the raw pattern.

## Marks

For each position: `cost = dG(wild type) - dG(deletion variant)`, within-domain statistics as
everywhere else in this project. Propensity scales fixed in advance: Pace & Scholtz 1998 helix
propensity (alanine 0.0, glycine 1.0, proline 3.16 - a penalty, so higher = worse helix former)
and Koehl-Levitt / Minor-Kim style beta-sheet propensity. Because I am not certain of the exact
published beta-sheet numbers from memory, the beta arm uses a **scale-free** form only.

- **D-1 (scale free, the primary mark).** An OLS of deletion cost on domain fixed effects, residue
  identity, structure class, relative solvent accessibility and residue volume is compared with
  the same model plus a residue-by-structure interaction. Mark: the F test for the interaction has
  p < 0.001, with standard errors clustered by domain.
- **D-2 (direction).** Within the interaction, the helix-versus-loop contrast correlates with Pace
  helix propensity across the 20 residues, Spearman <= -0.4 or >= 0.4 with the sign predicted
  below, and a bootstrap CI over domains excluding 0.
  Predicted sign: a residue that is a **poor** helix former (high Pace penalty) should be
  **cheaper** to delete from a helix than a good one, so the correlation between Pace penalty and
  (helix cost - loop cost) should be **negative**.
- **D-3 (not burial).** D-1 holds within each tertile of relative solvent accessibility.
- **D-4 (not volume).** D-1 holds after volume is included, and the helix contrast is not simply
  monotone in volume (report Spearman of the contrast against volume alongside propensity; if
  volume explains it better than propensity, H-del is not supported).
- **D-5 (replication).** The direction of D-2 holds in the 57 designed HEEH_KT miniproteins.

## Kill conditions

1. Permuting the deleted residue identity within domain must abolish D-1.
2. If fewer than 15 residue types have at least 100 deletions in both helices and loops, D-2 is
   reported as underpowered.
3. Glycine and proline are reported separately as well as included, since both are special cases
   that could drive a correlation on their own. If removing them abolishes D-2, that is stated.

## What I expect

D-1 is likely to pass simply because 20 residues times 3 classes is a lot of freedom in 50,000
measurements. D-2 is the real test, and I do not know its outcome. The honest risk is that any
interaction is driven by glycine and proline alone, which would make it a statement about two
special residues rather than about propensity, and kill condition 3 is there to catch that.
