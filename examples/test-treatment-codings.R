# ── Comparing Treatment Codings on CATE Recovery ─────────────────────────────
#
# This script fits CausalHorseForest under each treatment coding scheme
# (centered, binary, adaptive, invariant) on a simulated dataset with
# heterogeneous treatment effects and compares RMSE of the estimated CATE.
#
# Run from the package root:
#   Rscript tests/test-treatment-codings.R

library(ShrinkageTrees)

# ── 1. Simulate data ────────────────────────────────────────────────────────

set.seed(42)
n     <- 250
p     <- 5
X     <- matrix(rnorm(n * p), ncol = p)

# Propensity score depends on X1
e_x        <- pnorm(0.5 * X[, 1])
treatment  <- rbinom(n, 1, e_x)

# Prognostic function mu(x) = X1 + 0.5*X2
mu_x <- X[, 1] + 0.5 * X[, 2]

# Heterogeneous treatment effect: tau(x) = 1 + X3 + 0.5*X4
tau_x <- 1 + X[, 3] + 0.5 * X[, 4]

# Outcome: y = mu(x) + (Z - 0.5) * tau(x) + noise
sigma_true <- 1
y <- mu_x + (treatment - 0.5) * tau_x + rnorm(n, 0, sigma_true)

# ── 2. MCMC settings ────────────────────────────────────────────────────────

N_burn  <- 2000
N_post  <- 2000
n_trees <- 50

# ── 3. Fit each treatment coding ────────────────────────────────────────────

codings <- c("centered", "binary", "adaptive", "invariant")
results <- list()

for (coding in codings) {

  cat("Fitting treatment_coding =", coding, "... ")
  t0 <- proc.time()

  extra_args <- list()
  if (coding == "adaptive") {
    extra_args$propensity <- e_x
  }

  fit <- do.call(CausalHorseForest, c(list(
    y                        = y,
    X_train_control          = X,
    X_train_treat            = X,
    treatment_indicator_train = treatment,
    outcome_type             = "continuous",
    treatment_coding         = coding,
    number_of_trees          = n_trees,
    N_post                   = N_post,
    N_burn                   = N_burn,
    store_posterior_sample    = TRUE,
    verbose                  = FALSE
  ), extra_args))

  elapsed <- (proc.time() - t0)["elapsed"]

  # Estimated CATE = posterior mean of tau(x_i)
  tau_hat <- fit$train_predictions_treat

  # Posterior CATE samples (N_post x n matrix)
  tau_samples <- fit$train_predictions_sample_treat

  # RMSE of posterior mean CATE vs truth
  rmse <- sqrt(mean((tau_hat - tau_x)^2))

  # Bias of ATE
  ate_true <- mean(tau_x)
  ate_hat  <- mean(tau_hat)
  ate_bias <- ate_hat - ate_true

  # Posterior ATE samples
  ate_samples <- rowMeans(tau_samples)
  ate_ci      <- quantile(ate_samples, probs = c(0.025, 0.975))
  ate_covers  <- ate_true >= ate_ci[1] && ate_true <= ate_ci[2]

  # Acceptance ratio
  ar <- fit$acceptance_ratio_treat

  # Invariant coding: check b0, b1
  b0_mean <- NA; b1_mean <- NA
  if (coding == "invariant") {
    b0_mean <- mean(fit$b0)
    b1_mean <- mean(fit$b1)
  }

  results[[coding]] <- list(
    rmse       = rmse,
    ate_true   = ate_true,
    ate_hat    = ate_hat,
    ate_bias   = ate_bias,
    ate_ci     = ate_ci,
    ate_covers = ate_covers,
    ar_treat   = ar,
    b0_mean    = b0_mean,
    b1_mean    = b1_mean,
    elapsed    = elapsed
  )

  cat(sprintf("RMSE = %.3f, ATE bias = %+.3f, AR = %.3f (%.1fs)\n",
              rmse, ate_bias, ar, elapsed))
}

# ── 4. Summary table ────────────────────────────────────────────────────────

cat("\n")
cat("==============================================================\n")
cat("  Treatment Coding Comparison — CATE Recovery (n =", n, ")\n")
cat("==============================================================\n")
cat(sprintf("  True ATE = %.3f\n\n", mean(tau_x)))
cat(sprintf("  %-12s %8s %10s %8s %8s %6s\n",
            "Coding", "RMSE", "ATE bias", "AR", "Cover", "Time"))
cat("  ", strrep("-", 56), "\n", sep = "")

for (coding in codings) {
  r <- results[[coding]]
  cat(sprintf("  %-12s %8.3f %+10.3f %8.3f %8s %5.1fs\n",
              coding, r$rmse, r$ate_bias, r$ar_treat,
              ifelse(r$ate_covers, "yes", "no"), r$elapsed))
}

# Extra row for invariant b0, b1
r_inv <- results[["invariant"]]
if (!is.na(r_inv$b0_mean)) {
  cat(sprintf("\n  Invariant coding:  mean(b0) = %.4f,  mean(b1) = %.4f\n",
              r_inv$b0_mean, r_inv$b1_mean))
  cat("  (Expected: b0 ~ -0.5, b1 ~ 0.5 or similar non-degenerate values)\n")
}

cat("==============================================================\n")
