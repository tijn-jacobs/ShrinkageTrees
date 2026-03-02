# ShrinkageTrees x.x.x

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

## Other changes

- Fixed a latent bug in `HorseTrees` and `ShrinkageTrees` where `sigma_hat`, `y_mean`, and `lambda` were not initialised in the binary (probit) branch.

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
