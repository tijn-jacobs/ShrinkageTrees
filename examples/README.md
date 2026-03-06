# Examples

These scripts are for interactive testing and debugging during development.
They are **not** run by `devtools::test()` or `R CMD check` — those use the
automated tests in `tests/testthat/`.

**Note:** This directory is excluded from the built package tarball.
The scripts are available from the [GitHub repository](https://github.com/tijn-jacobs/ShrinkageTrees).

## Scripts

| Script | Description |
|---|---|
| `test-continuous.R` | HorseTrees / ShrinkageTrees on continuous outcomes, prior comparison, multi-chain |
| `test-binary.R` | HorseTrees on binary outcomes (probit BART), Brier score, AUC |
| `test-survival.R` | Right-censored and interval-censored survival, SurvivalBART wrapper, survival plots |
| `test-survival-curves.R` | Quick test of survival curve plotting and monotonicity (training + prediction) |
| `test-causal.R` | CausalHorseForest / CausalShrinkageForest, CATE RMSE, ATE recovery, multi-chain |
| `test-treatment-codings.R` | Comparison of all four treatment codings (centered, binary, adaptive, invariant) on CATE recovery |
| `test-plots.R` | Visual verification of all plot types, saves to PDF |
| `test-coda-diagnostics.R` | `as.mcmc.list()` S3 method, Gelman-Rubin, ESS, Geweke diagnostics |
| `test-ovarian.R` | Full worked example on the TCGA ovarian dataset: survival prediction + causal inference |

## Usage

Run any script from the package root:

```r
Rscript examples/test-continuous.R
```

Or source interactively in R:

```r
source("examples/test-causal.R")
```
