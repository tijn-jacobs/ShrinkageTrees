# ShrinkageTrees 2.1.0

## Two corrections to the horseshoe global update (breaking)

Two independent defects in `Horseshoe::GlobalUpdate()` have been fixed. Both
affected every fit using `prior_type = "horseshoe"` (including `HorseTrees()`
and `CausalHorseForest()`). Results from earlier versions will not reproduce
under this release.

1. **`global_hp` had no effect.** The scale entered the auxiliary draw as
   `1.0 / alpha_global * alpha_global`, which by operator precedence is the
   constant `1.0`. The user's `global_hp` was silently discarded, so every fit
   behaved as though `global_hp = 1` regardless of what was supplied. This was
   reported by a referee of the accompanying manuscript.

2. **The error scale and the forest scale were swapped.** `EtaPrior::GlobalUpdate()`
   declares `(..., sigma, omega, ...)` but was being called with `(omega, sigma)`.
   Both are `const double&`, so the swap compiled silently. The `tau` draw
   therefore divided by `sigma` where it should have divided by `omega`. The
   effect is small where `omega = 1` (`HorseTrees()`, `ShrinkageTrees()`) and
   substantial where `omega = 1/2` (the causal forests), where it grew with the
   number of covariates.

## Recalibrated default shrinkage (breaking)

The old default of `k = 0.1` was chosen while defect (1) was present, so it was
tuned for a model in which the global scale was pinned at 1. With the defect
corrected, that value shrinks far too aggressively: in our simulations it drove
pointwise coverage of the conditional treatment effect from roughly 0.95 down to
about 0.4. Re-running existing scripts against this release *without* updating
`k` is therefore worse than not upgrading at all.

The defaults are now calibrated for the corrected sampler:

| function | old | new |
|---|---|---|
| `HorseTrees()` | `k = 0.1` | `k = 1.0` |
| `CausalHorseForest()` | `k = 0.1` | `k = 1.5` |
| `ShrinkageTrees()` | `local_hp`, `global_hp` required | default to `1.0 / sqrt(number_of_trees)` |
| `CausalShrinkageForest()` | four scales required | default to `1.5 / sqrt(number_of_trees_*)` |

Values of `k` between roughly 0.5 and 1.5 (single-forest models) or 1 and 2
(causal models) worked well across a range of simulated settings, with smaller
values shrinking more aggressively and larger values being more conservative.
The defaults sit in the middle of each range. The
causal range is higher because the treatment forest enters as `b * tau(x)` with
`b = +/- 1/2`, halving its contribution to the response; measured on the scale
of that contribution the two recommendations nearly coincide. `k` remains
exposed on all four functions and is a natural target for cross-validation.

The `local_hp` / `global_hp` arguments are also no longer mandatory for
`prior_type = "horseshoe"`; either may be supplied alone, with the other taking
its default. Only their product identifies the prior, so setting them equal
loses nothing. `prior_type = "horseshoe_fw"` is unchanged and still requires
both to be given explicitly.

## `predict()` passed the training design in the wrong memory order

`predict()` on a `ShrinkageTrees` object handed the stored training design
matrix straight to the C++ sampler. The C++ layer reads designs as a flat
**row-major** buffer, but R coerces a matrix to a vector **column-major**, so
the training data was transposed inside the sampler on every call whenever
`n != p`. The test design on the same call was flattened correctly, so the two
disagreed. It ran and returned plausible numbers.

`predict()` on a `CausalShrinkageForest` object was unaffected: it already
flattened both designs correctly.

## New: `posterior_projection()`

Lower-dimensional posterior summarisation following Woody, Carvalho and Murray
(2021). Every posterior draw of the fitted function, or of the treatment-effect
and prognostic surfaces for causal fits, is projected onto a simpler summary
model: a linear model (unpenalised, ridge, lasso, or elastic net), a
spline-additive model, or a shallow CART tree. Returns a posterior over
summaries with credible intervals, plus the posterior of the summary R-squared.
`glmnet` and `rpart` are new suggested packages.

## Reproducing earlier output

There is no setting that reproduces pre-2.1.0 results exactly, because the
sampler itself has changed. The closest correspondence is that the old code was
equivalent to the corrected code with `global_hp = 1`, since the global scale
was inert.

# ShrinkageTrees 2.0.2

## Bayesian bootstrap for the average treatment effect

`summary()`, `plot(type = "ate")`, and `predict()` for
`CausalShrinkageForest` (and `CausalHorseForest`) now default to a
Bayesian-bootstrap posterior for the average treatment effect: at each
MCMC iteration the per-observation CATEs are reweighted with
Dirichlet(1, ..., 1) weights before being summed, giving a draw from the
posterior of the *population* ATE (PATE). Credible intervals are
correspondingly wider than before because they now propagate uncertainty
in the covariate distribution, not only in `tau(x)`.

- New argument `bayesian_bootstrap = TRUE` on `summary.CausalShrinkageForest()`,
  `plot.CausalShrinkageForest()` (for `type = "ate"`), and
  `predict.CausalShrinkageForest()`. Set to `FALSE` to recover the
  previous equal-weight mixed ATE (MATE).
- `predict()` now also returns an `ate` summary for the new data and
  retains the posterior CATE sample matrix as `cate_samples`.
- New exported helper `bayesian_bootstrap_ate()` returns PATE and MATE
  summaries (means, CIs, and full posterior draws) from either a fitted
  `CausalShrinkageForest` or a `CausalShrinkageForestPrediction`.

This is a *breaking change for printed/plotted numerics*: existing
scripts will report wider CIs than before. Use
`bayesian_bootstrap = FALSE` to reproduce previous output.

# ShrinkageTrees 2.0.1

## `ovarian` dataset restructured

The `ovarian` dataset is now a single data frame (previously a list with
`$clinical` and `$X` elements). Clinical columns (`OS_time`, `OS_event`,
`treatment`, `age`, `figo_stage`, `tumor_grade`) and the 2000 gene
expression columns are combined into one data frame with 2006 columns.
This simplifies data access and aligns the format with the `pdac` dataset.

Code that previously used `ovarian$clinical` or `ovarian$X` must be
updated — see `?ovarian` for the new structure.

## Bug fixes

- Fixed `plot(fit, type = "ate")` and `plot(fit, type = "cate")` for
  `CausalShrinkageForest` models incorrectly subtracting the control forest
  predictions from the treatment forest predictions.

# ShrinkageTrees 2.0.0

## TCGA ovarian cancer dataset (`ovarian`)

Added the `ovarian` dataset: a processed TCGA-OV cohort (n = 357) for
high-dimensional survival prediction and causal inference.

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
