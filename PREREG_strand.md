# Pre-registration 3: are beta strands more indel-tolerant than alpha helices within a protein?

Written 26 Nov 2026, after the decomposition in 28/29 and **before** running the controls and the
replication below.

## The observation

Within the 84 Tsuboyama domains that contain both helix and strand, paired within each domain:

| | strand minus helix | 95% CI by domain |
|---|---|---|
| single-Ala insertion cost | **-0.478** | [-0.665, -0.287] |
| single-residue deletion cost | **-0.445** | [-0.669, -0.227] |

Negative = cheaper in the strand. Relative to each domain's own loops, insertion costs +0.359
[0.238, 0.483] in strands against +0.837 [0.677, 0.996] in helices.

Pooled over all domains without controlling for domain, the sign reverses: insertion costs 2.191
in strands against 1.707 in helices, because beta-rich domains are globally more indel-sensitive
(mean insertion cost 2.444 in all-beta domains, 1.386 in all-alpha). The published claims point
the other way: Yang et al. 2024 report strands as the most indel-depleted class in evolution (O/E
0.43 against 0.57 for helices), and Topolska et al. 2025 report strands as slightly MORE sensitive
than helices. If their comparisons are made across proteins, this reversal is the explanation.

**H-strand: within a protein, a single-residue indel costs less in a beta strand than in an alpha
helix.**

## Why this could still be wrong, and what must be shown

Each of these is a mandatory control, not an optional one. H-strand is only reported as supported
if **all five** hold.

- **C1 burial.** Strands and helices differ in solvent exposure. The paired difference must keep
  its sign within every tertile of relative solvent accessibility.
- **C2 confidence.** Same, within every tertile of AlphaFold pLDDT.
- **C3 element position.** My helix definition requires a run of 5 with two consecutive 4-turns,
  the strand definition only a 2-residue ladder, so helix positions may sit deeper inside their
  element. The difference must survive comparing positions matched on distance to the element
  terminus (middles against middles, ends against ends).
- **C4 residue composition.** Strands are rich in beta-branched and hydrophobic residues. The
  difference must survive restricting to positions with the same wild-type residue, and must hold
  for a majority of individual residue types.
- **C5 independent structure labels.** The difference must hold when the Lehner lab's own
  secondary-structure table is used instead of 17_dssp.py.

## Replication, pre-registered

Topolska et al.'s yeast abundancePCA data (9 domains, 1,420 insertions, 1,450 deletions, their own
structure labels, orthogonal assay).

- **RS1.** The within-domain paired strand-minus-helix insertion cost is negative, in at least
  **6 of** the domains that contain both classes, with a by-domain CI excluding 0.
- **RS2.** Same for deletions.

Note in advance: only 9 domains exist, several may lack one class, and the earlier asymmetry
replication (PREREG_bulge.md R1) **failed** on this dataset. Power is low. If RS1 fails while C1-C5
pass, the honest conclusion is a strong single-dataset result awaiting replication, and it will be
written that way.

## Added 26 Nov 2026, before running: a second replication set found inside Lib3

Lib3 also contains single indels in **designed** miniproteins that neither Topolska et al. nor
Gutierrez & Rocklin appear to have used: 69 HEEH_KT designs, 57 of which have an AlphaFold model
here, with about 131 single-indel variants each. Their topology is helix-strand-strand-helix by
construction, so every one of them supplies a within-protein helix-versus-strand contrast, and
they share no evolutionary history with the natural domains that produced the observation. The
numbering will be verified by checking that the residue named in each deletion variant matches
the model at that position; if the match rate is below 0.95 the set is dropped.

- **RS3.** The within-protein paired strand-minus-helix insertion cost is negative, with a
  by-domain CI excluding 0 over at least 30 designs.
- **RS4.** Same for deletions.

These marks are set before looking. RS3/RS4 are the strongest available test of H-strand, because
the designs are numerous, structurally uniform and independent of the discovery set.

## Kill conditions

1. Permutation of cost within domain must abolish the paired difference.
2. If the difference is carried by fewer than 30 domains, it is not reported.
3. If C3 shows the effect only at element ends, it is a boundary artefact of the structure
   assignment and H-strand is withdrawn.

## What I expect

The effect is large and the within-domain pairing already removes the obvious confounder, so I
give H-strand better than even odds on C1-C5. I expect RS1 to be underpowered rather than
contradictory. The claim that would be genuinely new is not "strands tolerate indels" on its own
but the **reversal**: that the pooled comparison in the literature has the sign backwards because
it mixes proteins.
