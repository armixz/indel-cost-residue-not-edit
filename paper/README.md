# Paper build

```
# 1. analyses (WSL, /root/envs/indel/bin/python)
python ../analysis/41_verify.py  <lib3.csv> <afdir> ../results ../results/verify.json
python ../analysis/42_predictive.py <lib3.csv> <afdir> ../results ../results/predictive.json
python ../analysis/43_discriminate.py <lib3.csv> <afdir> ../results ../results/discriminate.json
python ../analysis/44_which_term.py <lib3.csv> <afdir> ../results/which_term.json

# 2. figures (vector PDF)
python make_figures.py <lib3.csv> <afdir> ../results .

# 3. numbers -> numbers.tex  (fails loudly if any macro has no source)
python update_numbers.py ../results . numbers.tex

# 4. manuscript
latexmk -pdf main.tex
```

`numbers.tex` is generated, never edited by hand: every number quoted in `main.tex` is a macro
defined there and sourced from an analysis output, so the text cannot drift from the analysis.
A macro with no source compiles as a red **MISSING**, so an incomplete build is visible.

Pre-registrations for every test, including the ones whose marks failed and the two framings
withdrawn during the work, are in `../PREREG_*.md`. The running record with its own correction
section is `../RESULTS_indel.md`.
