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
| `test-bayesian-bootstrap-ate.R` | Bayesian-bootstrap PATE against the plug-in MATE |
| `semi-synthesise-ovarian.R` | Rebuilds the shipped `ovarian` and `ovarian_truth` datasets, with diagnostics |

`semi-synthesise-ovarian.R` is not a test. It is the generating script for the
two shipped datasets: it reads the installed `ovarian` for its real covariates,
resimulates treatment, the confounder and the outcomes from a known process,
plots the diagnostics, and writes `data/*.rda` only when called with `--save`.
It fits no models; the analysis of this data is
`simulations/r-journal-paper/ovarian_analysis.R`.

```r
Rscript examples/semi-synthesise-ovarian.R           # diagnostics only
Rscript examples/semi-synthesise-ovarian.R --save    # overwrite data/*.rda
```

## Usage

Run any script from the package root:

```r
Rscript examples/test-continuous.R
```

Or source interactively in R:

```r
source("examples/test-causal.R")
```

