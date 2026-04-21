# Bayesian Bootstrap for PATE: Implementation Note

## Overview

Add a post-hoc function that takes the existing posterior CATE samples from a fitted `CausalShrinkageForest` object and reweights them with Bayesian bootstrap weights to produce PATE estimates with properly calibrated credible intervals.

This is purely post-hoc: no changes to the MCMC sampler are needed. The posterior CATE samples $\tau^{(s)}(x_i)$ are already stored; we just change how we aggregate them into an ATE.

## What we have

After fitting a causal model with `store_posterior_sample = TRUE`, the fitted object contains a matrix of posterior CATE draws:

- Rows: MCMC iterations $s = 1, \ldots, S$
- Columns: observations $i = 1, \ldots, n$

Currently the ATE is computed as:

$$\widehat{\text{ATE}}^{(s)} = \frac{1}{n} \sum_{i=1}^{n} \tau^{(s)}(x_i)$$

This is the MATE. It conditions on the observed covariates as fixed.

## What we want

At each posterior draw $s$, draw Dirichlet weights and compute:

$$\widehat{\text{PATE}}^{(s)} = \sum_{i=1}^{n} w_i^{(s)} \, \tau^{(s)}(x_i), \quad (w_1^{(s)}, \ldots, w_n^{(s)}) \sim \text{Dir}(1, \ldots, 1)$$

The collection $\{\widehat{\text{PATE}}^{(s)}\}_{s=1}^{S}$ gives the posterior distribution of the PATE on the log-survival scale. Credible intervals from this distribution incorporate uncertainty in both $\tau(\cdot)$ and $F_X$.

## Implementation plan

### Function signature

```r
bayesian_bootstrap_ate <- function(object, alpha = 0.05, seed = NULL)
```

- `object`: a fitted `CausalShrinkageForest` with `store_posterior_sample = TRUE`
- `alpha`: credible interval level (default 0.05 for 95% CI)
- `seed`: optional seed for reproducibility of the Dirichlet draws
- Returns: a list with `pate_mean`, `pate_ci`, `pate_samples`, `mate_mean`, `mate_ci`

### Core logic (R, ~15 lines)

```r
bayesian_bootstrap_ate <- function(object, alpha = 0.05, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  # tau_samples: S x n matrix of posterior CATE draws
  tau_samples <- object$tau_posterior_samples  # or however it's stored
  S <- nrow(tau_samples)
  n <- ncol(tau_samples)
  
  # MATE (current approach, for comparison)
  mate_samples <- rowMeans(tau_samples)
  
  # PATE via Bayesian bootstrap
  # Draw S sets of Dirichlet(1,...,1) weights
  # Dirichlet(1,...,1) = normalised Exp(1) draws
  W <- matrix(rexp(S * n), nrow = S, ncol = n)
  W <- W / rowSums(W)
  
  # Weighted ATE at each iteration
  pate_samples <- rowSums(W * tau_samples)
  
  # Summaries
  list(
    pate_mean = mean(pate_samples),
    pate_ci   = quantile(pate_samples, probs = c(alpha / 2, 1 - alpha / 2)),
    pate_samples = pate_samples,
    mate_mean = mean(mate_samples),
    mate_ci   = quantile(mate_samples, probs = c(alpha / 2, 1 - alpha / 2)),
    mate_samples = mate_samples
  )
}
```

### Key details

- **Dirichlet sampling**: `rexp(S * n)` draws S × n independent Exp(1) values. Normalising each row to sum to 1 gives Dir(1,...,1). This is the standard trick; no extra packages needed.
- **Scale**: Everything stays on the log-survival scale. The output is the PATE of $\tau(x) = E[\log T(1) - \log T(0) \mid X = x]$ integrated over $F_X$. No back-transformation.
- **Memory**: The `W` matrix is S × n (e.g. 20000 × 357 ≈ 57MB for the ovarian example). This is fine. If memory is a concern, loop over iterations instead of materialising the full matrix.

### What needs checking in the package

1. **Where are the CATE samples stored?** Find the field name in the `CausalShrinkageForest` object that holds the S × n matrix of posterior $\tau(x_i)$ draws. It may be `object$tau_posterior_samples`, `object$cate_samples`, or similar. If `store_posterior_sample = FALSE` (the default), the samples are discarded and only the posterior mean is kept — in that case the function should error with a clear message.

2. **Are samples pooled across chains?** If `n_chains > 1`, the samples from all chains should already be concatenated row-wise (the paper says they are). Verify this.

3. **Where to put the function.** Add as an exported function in the package, or as an S3 method like `bayesian_bootstrap.CausalShrinkageForest`. A standalone function is simpler.

### Testing

- Verify that `pate_ci` is wider than `mate_ci` (it should always be).
- Verify that `pate_mean ≈ mate_mean` (they should be close; the Dirichlet weights have mean 1/n).
- With `seed` set, output should be reproducible.

### For the paper

Once implemented, rerun the ovarian causal analysis and report both:

```r
bb <- bayesian_bootstrap_ate(fit_causal, seed = 42)
cat("MATE:", round(bb$mate_mean, 4), "95% CI:", round(bb$mate_ci, 4), "\n")
cat("PATE:", round(bb$pate_mean, 4), "95% CI:", round(bb$pate_ci, 4), "\n")
```

This directly addresses the known limitation in the Discussion section and turns it from "future work" into a solved problem.
