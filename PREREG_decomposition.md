# Pre-registration 4: is a deletion an insertion plus the loss of a side chain?

Written 26 Nov 2026, before extracting any substitution data from Lib3.

## The question

Fitting deletion cost against insertion cost within each protein gives slope 0.692 [0.664, 0.720]
and intercept **0.81 [0.71, 0.92] kcal/mol**: even where an insertion is free, a deletion is not.
What is that offset?

A deletion does two things at once - it removes the residue's side chain and it perturbs the
backbone. An insertion does only the second. A substitution to glycine does only the first. The
three measurements exist at the same positions in the same library, so the decomposition is
directly testable.

**H-decomp: cost(deletion at i) = cost(insertion at i) + cost(X_i -> Gly at i), up to a constant.**

If it holds, the insertion/deletion asymmetry is not a separate phenomenon at all: it is the
side-chain term, and everything found so far follows from it. In particular:

- the asymmetry should be large where the side chain contributes a lot (buried, large residues)
  and near zero where it does not (glycine positions, exposed loops)
- the beta-strand excess should follow from strands having side chains that contribute more
- the headroom pattern should follow from the two terms saturating together

If it fails, the asymmetry is genuinely about the backbone, and the bulge-style explanations stay
on the table.

## Marks

Let `d` = deletion cost, `i` = insertion cost (Ala), `g` = cost of substituting the wild-type
residue to glycine, all at the same position, all within the same protein.

- **D1 (the core).** Within domain, regressing `d` on `i + g` gives a slope whose CI contains 1.0,
  and the residual RMSE is at least 15% lower than regressing `d` on `i` alone.
- **D2 (the offset is absorbed).** The intercept of `d ~ i + g` is smaller in magnitude than the
  0.81 kcal/mol found for `d ~ i`, with a by-domain CI that excludes 0.81.
- **D3 (glycine positions).** At positions whose wild-type residue is already glycine, `g` is zero
  by construction, so the asymmetry `d - i` should be near zero. Mark: mean `d - i` at
  wild-type-glycine positions is below 0.2 kcal/mol in magnitude, against 0.303 [0.257, 0.351]
  overall.
- **D4 (the strand excess).** After subtracting `g`, the strand-minus-loop asymmetry should shrink
  substantially. Mark: it falls below half of -0.324, i.e. above -0.162. **If instead it survives
  at close to full size, the strand effect is a backbone effect and H-decomp does not explain it.**

D3 and D4 are the informative ones: D3 is a near-parameter-free prediction, and D4 decides whether
the strand result is side-chain or backbone in origin.

## Controls

1. Use the same quality filter (dG CI <= 1 kcal/mol) for substitutions as for indels.
2. Replicate every mark on the designed HEEH_KT miniproteins.
3. Permutation: shuffling `g` within domain must break D1.
4. Report how many positions have all three measurements; if fewer than 3,000 across fewer than
   100 domains, the test is reported as underpowered.

## What I expect

D3 seems likely to pass at least partly - glycine has no side chain to lose. D4 is the one I
genuinely do not know. If the strand effect is bulge-related it is a backbone effect and should
survive; if it is side-chain packing in a sheet it should vanish. Either answer is informative and
will be reported as it falls.
