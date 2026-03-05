# Manual Tests

These scripts are for interactive testing and debugging during development.
They are **not** run by `devtools::test()` or `R CMD check` — those use the
automated tests in `tests/testthat/`.

## Scripts

| Script | Description |
|---|---|
| `test-continuous.R` | HorseTrees / ShrinkageTrees on continuous outcomes, prior comparison, multi-chain |
| `test-binary.R` | HorseTrees on binary outcomes (probit BART), Brier score, AUC |
| `test-survival.R` | Right-censored and interval-censored survival, SurvivalBART wrapper, survival plots |
| `test-causal.R` | CausalHorseForest / CausalShrinkageForest, CATE RMSE, ATE recovery, multi-chain |
| `test-treatment-codings.R` | Comparison of all four treatment codings (centered, binary, adaptive, invariant) on CATE recovery |
| `test-plots.R` | Visual verification of all plot types, saves to PDF |

## Usage

Run any script from the package root:

```r
Rscript tests/manual/test-continuous.R
```

Or source interactively in R:

```r
source("tests/manual/test-causal.R")
```
