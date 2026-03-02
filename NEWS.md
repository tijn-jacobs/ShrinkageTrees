# ShrinkageTrees 1.3.0

## Multi-chain MCMC (`n_chains`)

All four primary model-fitting functions — `ShrinkageTrees`, `HorseTrees`,
`CausalHorseForest`, and `CausalShrinkageForest` — now accept an `n_chains`
argument (default `1`). When `n_chains > 1`:

- Independent chains are dispatched in parallel using `parallel::mclapply`
  (falls back to sequential execution on Windows).
- The number of cores used is `min(n_chains, parallel::detectCores())`.
- Posterior sample matrices from all chains are row-bound into a single matrix,
  giving `N_post * n_chains` total draws.
- Posterior means are recomputed from the pooled samples; sigma and acceptance
  ratio vectors are concatenated across chains.
- The returned object is a standard `ShrinkageTrees` or `CausalShrinkageForest`
  instance, so all existing `print`, `summary`, and `predict` methods work
  without modification.
- The survival wrappers `SurvivalBART`, `SurvivalDART`, `SurvivalBCF`, and
  `SurvivalShrinkageBCF` inherit `n_chains` support through `...`.
- `print` and `summary` output adapts automatically: single-chain models show
  _Posterior draws_, multi-chain models show _Chains_ and _Draws per chain_,
  with per-chain acceptance ratios listed separately.

## S3 classes and methods

- Added S3 classes `ShrinkageTrees` and `CausalShrinkageForest` with constructors in `constructors.R`.
- Added `print` methods for both classes, displaying model specification, MCMC settings, acceptance ratio, and posterior mean sigma.
- Added `summary` methods for both classes, returning an inspectable object with posterior sigma (mean, SD, 95% CI), prediction summaries, variable importance (posterior inclusion probabilities), and — for causal models — ATE with credible interval (when `store_posterior_sample = TRUE`) and CATE heterogeneity.
- Added `predict` method for `ShrinkageTrees`, enabling posterior predictive inference on new data by re-running the sampler with stored training data and hyperparameters. Returns a `ShrinkageTreesPrediction` object with posterior mean and credible interval bounds.
- Added `print` and `summary` methods for `ShrinkageTreesPrediction`, showing a formatted head-style table and a min/quartile/max distribution summary respectively.

## Posterior visualisation (`plot`)

S3 `plot()` methods added for `ShrinkageTrees` and `CausalShrinkageForest`.
Requires the suggested packages `bayesplot` and `ggplot2`.

- `plot(fit, type = "trace")` — sigma traceplot; one line per chain, useful for assessing mixing.
- `plot(fit, type = "density")` — overlaid posterior density of sigma, one curve per chain.
- `plot(fit, type = "vi")` — posterior credible intervals for variable inclusion probabilities (top `n_vi` predictors).
- `plot(fit, type = "ate")` — posterior density of the ATE with 95 % credible region _(causal models only; requires `store_posterior_sample = TRUE`)_.
- `plot(fit, type = "cate")` — point estimates and 95 % credible intervals for the CATE of each training observation, sorted by posterior mean _(causal models only; requires `store_posterior_sample = TRUE`)_.
- `plot(fit, type = "vi", forest = "both")` — side-by-side VI for the control and treatment forests _(causal models only)_.

## Vignette

- Added a package vignette (_ShrinkageTrees: Introduction and Usage_)
  demonstrating all main functions (`HorseTrees`, `ShrinkageTrees`,
  `SurvivalBART`, `SurvivalDART`, `SurvivalBCF`, `SurvivalShrinkageBCF`,
  `CausalHorseForest`, `CausalShrinkageForest`), all S3 methods (`print`,
  `summary`, `predict`, `plot`), multi-chain MCMC, and a full TCGA PAAD
  case study.

## Bug fixes

- Fixed `CausalHorseForest` and `CausalShrinkageForest` failing with
  `"argument 'y_train' is missing"` when called directly (broken constructor
  call introduced in 1.3.0 S3 refactor).
- Fixed `plot(..., type = "vi")` crashing with `"argument must be coercible
  to non-negative integer"` in all four model functions: covariate matrices
  were being stored as flat numeric vectors instead of matrices, making
  `ncol()` return `NULL`.
- Fixed a latent bug in `HorseTrees` and `ShrinkageTrees` where `sigma_hat`,
  `y_mean`, and `lambda` were not initialised in the binary (probit) branch.

# ShrinkageTrees 1.2.0

- Added `SurvivalBCF` wrapper for AFT-based Bayesian Causal Forests.
- Added `SurvivalShrinkageBCF` with Dirichlet structural sparsity.
- Added `SurvivalDART` for sparse high-dimensional survival modeling.
- Added `SurvivalBART` for survival modeling using standard BART.
- Unified survival wrappers with flexible argument forwarding via `...`.

# ShrinkageTrees 1.1.0

- Replaced the internal `std::map` structure with a more efficient vector-based lookup, improving overall computational speed by approximately 30%.
- Added the `standard` prior type option for `ShrinkageTrees`, corresponding to the conventional BART implementation without reversible jumps.
- Added missing checks for training data dimensions to improve input validation and error handling.

# ShrinkageTrees 1.0.3

- Refactored the non-reversible tree modification routines in the C++ backend for improved clarity and maintainability.
- Corrected the 'leafs' typo throughout the codebase (now consistently 'leaves').

# ShrinkageTrees 1.0.2

- Improved handling of censored survival outcomes in the back-end
  (preparatory changes for interval censoring support).
- Minor internal refactoring; no changes to the user-facing API.

# ShrinkageTrees 1.0.1

- Fixed bugs in the demo script
- Added a `CONTRIBUTING.md` file with guidelines for contributors
- Corrected minor typos in the source code (non-functional changes)

# ShrinkageTrees 1.0.0

🎉 First CRAN release of **ShrinkageTrees**!

This package provides Bayesian regression tree models with shrinkage priors, supporting:

- Continuous outcomes
- Binary outcomes
- Right-censored survival data

It includes four core functions:

- `HorseTrees()`: fits a single regression tree with a standard Horseshoe prior.
- `ShrinkageTrees()`: fits a single tree with customizable shrinkage priors.
- `CausalHorseForest()`: fits a causal forest using the standard Horseshoe prior.
- `CausalShrinkageForest()`: fits a flexible causal forest with user-defined shrinkage priors and tuning options.

The `...Trees` functions use a single learner to estimate the outcome model directly. In contrast, the `Causal...Forest` variants fit separate models for the treated and control regression function. This enables estimation of conditional average treatment effects (CATEs).
