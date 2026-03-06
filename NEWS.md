# ShrinkageTrees 2.0.0

## TCGA ovarian cancer dataset (`ovarian`)

Added the `ovarian` dataset: a processed TCGA-OV cohort (n = 357) for
high-dimensional survival prediction and causal inference. The dataset is
a list with two elements:

- `X`: a 357 x 2000 gene expression matrix (log2-normalised TPM, top 2000
  most variable genes by median absolute deviation).
- `clinical`: a data frame with overall survival time and event indicator,
  age, FIGO stage, tumor grade, and a binary treatment indicator
  (carboplatin vs cisplatin).

See `?ovarian` and `examples/test-ovarian.R` for a full worked example
covering survival prediction (SurvivalBART, SurvivalDART, HorseTrees) and
causal inference (SurvivalBCF, SurvivalShrinkageBCF, CausalHorseForest).

## Treatment coding for causal models (`treatment_coding`)

All causal model functions — `CausalHorseForest()`, `CausalShrinkageForest()`,
`SurvivalBCF()`, and `SurvivalShrinkageBCF()` — now accept a `treatment_coding`
argument controlling how the treatment indicator enters the BCF decomposition
y = f(x) + b \* tau(x) + epsilon. Four options are available:

- `"centered"` (default): b_i in {-1/2, 1/2}. This is the original behaviour.
- `"binary"`: b_i in {0, 1}. Standard binary coding.
- `"adaptive"`: b_i = z_i - e_hat(x_i), where e_hat(x_i) is the estimated
  propensity score. This follows Hahn, Murray & Carvalho (2020) and is
  implemented in the `bcf` package. Requires a `propensity` vector.
- `"invariant"`: Parameter-expanded (invariant) treatment coding. The coding
  parameters b_0 and b_1 are assigned N(0, 1/2) priors and estimated within
  the Gibbs sampler via conjugate normal updates, yielding a parameterisation
  that is invariant to the coding of the treatment indicator (Hahn et al.,
  2020, Section 5.2). The treatment effect is tau(x) = (b_1 - b_0) \* tau_tilde(x).
  Posterior draws of b_0 and b_1 are returned in the fitted object.

The `predict()` method for `CausalShrinkageForest` objects automatically
carries forward the treatment coding used at training time. A `propensity_test`
argument is available for supplying test-set propensity scores (defaults to 0.5).

## Interval-censored survival outcomes

All survival-capable functions now support **interval-censored** data in
addition to right-censored data. Supply `left_time` and `right_time`
vectors (with `outcome_type = "interval-censored"`) instead of `y` and
`status`. Three censoring types are distinguished:

- **Exact events**: `left_time == right_time`.
- **Interval-censored**: finite `left_time < right_time`.
- **Right-censored**: `right_time = Inf`.

This convention follows `survival::Surv(type = "interval2")`. Censored
event times are augmented within the
AFT Gibbs sampler. The following functions are affected:

- `HorseTrees()`, `ShrinkageTrees()` (single-forest models)
- `CausalHorseForest()`, `CausalShrinkageForest()` (causal models)
- `SurvivalBART()`, `SurvivalDART()`, `SurvivalBCF()`,
  `SurvivalShrinkageBCF()` (survival wrappers)

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
- Added `predict` method for `ShrinkageTrees`, enabling posterior predictive inference on new data by re-running the sampler with stored training data and hyperparameters. Returns a `ShrinkageTreesPrediction` object with posterior mean and credible interval bounds. For survival models, the prediction object additionally stores `predictions_sample` (full posterior draws on the original scale) and `sigma` (posterior draws on the log-time scale), enabling posterior predictive survival curve plotting.
- Added `predict` method for `CausalShrinkageForest`, returning a `CausalShrinkageForestPrediction` object with three components: `prognostic` ($\mu(X)$), `cate` ($\tau(X)$), and `total` ($\mu(X) + \tau(X)$), each with posterior mean and credible interval bounds. For survival models with `timescale = "time"`, predictions are back-transformed to the original time scale and the CATE becomes a multiplicative time ratio.
- Added `print` and `summary` methods for `ShrinkageTreesPrediction` and `CausalShrinkageForestPrediction`.

## MCMC convergence diagnostics (`coda`)

- Added `as.mcmc.list()` S3 method for `ShrinkageTrees` objects, converting
  the sigma posterior (split by chain) into a `coda::mcmc.list`. This enables
  all standard `coda` diagnostics: Gelman–Rubin R-hat, effective sample size,
  Geweke test, Heidelberger–Welch test, autocorrelation plots, and more.
- `summary()` now automatically reports **effective sample size** (ESS) and —
  for multi-chain fits — the **Gelman–Rubin R-hat** when the suggested package
  `coda` is installed.
- Added `coda` to `Suggests` in DESCRIPTION.

## Posterior visualisation (`plot`)

S3 `plot()` methods added for `ShrinkageTrees`, `CausalShrinkageForest`,
and `ShrinkageTreesPrediction`. Requires the suggested package `ggplot2`.

- `plot(fit, type = "trace")` — sigma traceplot; one line per chain, useful for assessing mixing.
- `plot(fit, type = "density")` — overlaid posterior density of sigma, one curve per chain.
- `plot(fit, type = "vi")` — posterior credible intervals for variable inclusion probabilities (top `n_vi` predictors).
- `plot(fit, type = "ate")` — posterior density of the ATE with 95 % credible region _(causal models only; requires `store_posterior_sample = TRUE`)_.
- `plot(fit, type = "cate")` — point estimates and 95 % credible intervals for the CATE of each training observation, sorted by posterior mean _(causal models only; requires `store_posterior_sample = TRUE`)_.
- `plot(fit, type = "vi", forest = "both")` — side-by-side VI for the control and treatment forests _(causal models only)_.
- `plot(fit, type = "survival")` — posterior survival curves
  $S(t | x_i) = 1 - \Phi((\log t - \mu_i) / \sigma)$ derived from the AFT
  log-normal model _(survival outcomes only)_.
  - **Population-averaged curve** (default, `obs = NULL`): computes
    $\bar{S}(t) = n^{-1} \sum_i S(t | x_i)$ at each MCMC iteration with
    pointwise credible bands.
  - **Individual curves** (`obs = c(1, 5, ...)`): one curve per selected
    training observation with its own credible band.
  - `level` controls the credible band width (default 0.95).
  - `t_grid` allows a custom time grid; auto-generated if `NULL`.
  - `km = TRUE` overlays the Kaplan–Meier estimate as a dashed black
    step function (population-averaged plot only; requires `survival`
    package). Ignored with a message when `obs` is not `NULL`.
- `plot(pred, type = "survival")` — posterior **predictive** survival curves
  for new (out-of-sample) data from `predict()`. Same `obs`, `t_grid`, and
  `level` arguments as above. The KM overlay is not available for prediction
  objects.

## Vignette

- Added a package vignette (_ShrinkageTrees: Introduction and Usage_)
  demonstrating all main functions (`HorseTrees`, `ShrinkageTrees`,
  `SurvivalBART`, `SurvivalDART`, `SurvivalBCF`, `SurvivalShrinkageBCF`,
  `CausalHorseForest`, `CausalShrinkageForest`), all S3 methods (`print`,
  `summary`, `predict`, `plot`), multi-chain MCMC, and a full TCGA PAAD
  case study.

## Survival wrapper improvements

- `SurvivalBART()`, `SurvivalDART()`, `SurvivalBCF()`, and
  `SurvivalShrinkageBCF()` now accept `store_posterior_sample` as an
  explicit parameter (default `TRUE`), avoiding a "matched by multiple
  actual arguments" error when passing it via `...`.

## Bug fixes

- Fixed `CausalHorseForest` and `CausalShrinkageForest` failing with
  `"argument 'y_train' is missing"` when called directly (broken constructor
  call introduced in 2.0.0 S3 refactor).
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
