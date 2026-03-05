# ── Manual test: causal forest models ─────────────────────────────────────────
#
# Tests CausalHorseForest and CausalShrinkageForest on simulated data
# with known heterogeneous treatment effects. Checks ATE recovery, CATE
# RMSE, posterior samples, predict(), and plot().
#
# Run from the package root:
#   Rscript tests/manual/test-causal.R

library(ShrinkageTrees)

# ── 1. Data generation ──────────────────────────────────────────────────────

set.seed(42)
n <- 200; p <- 5
X <- matrix(rnorm(n * p), ncol = p)
X_test <- matrix(rnorm(50 * p), ncol = p)

treatment      <- rbinom(n, 1, 0.5)
treatment_test <- rbinom(50, 1, 0.5)

# Prognostic function
mu_x <- X[, 1] + 0.5 * X[, 2]

# Heterogeneous treatment effect
tau_x      <- 1 + X[, 3] + 0.5 * X[, 4]
tau_x_test <- 1 + X_test[, 3] + 0.5 * X_test[, 4]

# Outcome (centered coding: b = Z - 0.5)
y <- mu_x + (treatment - 0.5) * tau_x + rnorm(n)

cat(sprintf("True ATE: %.3f\n", mean(tau_x)))

# ── 2. CausalHorseForest ────────────────────────────────────────────────────

cat("\n=== CausalHorseForest ===\n")

set.seed(1)
fit_chf <- CausalHorseForest(
  y = y,
  X_train_control = X,
  X_train_treat = X,
  treatment_indicator_train = treatment,
  X_test_control = X_test,
  X_test_treat = X_test,
  treatment_indicator_test = treatment_test,
  outcome_type = "continuous",
  number_of_trees = 50,
  N_post = 1000, N_burn = 500,
  store_posterior_sample = TRUE,
  verbose = TRUE
)

cat("\n")
print(fit_chf)
cat("\n")
smry <- summary(fit_chf)

# CATE RMSE
tau_hat <- fit_chf$train_predictions_treat
rmse_cate <- sqrt(mean((tau_hat - tau_x)^2))
cat(sprintf("CATE RMSE (train): %.3f\n", rmse_cate))

# ATE
ate_hat <- mean(tau_hat)
cat(sprintf("ATE estimate: %.3f (true: %.3f)\n", ate_hat, mean(tau_x)))

# Acceptance ratios
cat(sprintf("AR control: %.3f\n", fit_chf$acceptance_ratio_control))
cat(sprintf("AR treat:   %.3f\n", fit_chf$acceptance_ratio_treat))

# Predict on new data
pred <- predict(fit_chf, newdata_control = X_test, newdata_treat = X_test,
                treatment_indicator = treatment_test)
cat(sprintf("Predict test CATE RMSE: %.3f\n",
            sqrt(mean((pred$mean_treat - tau_x_test)^2))))

# ── 3. CausalShrinkageForest ───────────────────────────────────────────────

cat("\n=== CausalShrinkageForest ===\n")

set.seed(1)
fit_csf <- CausalShrinkageForest(
  y = y,
  X_train_control = X,
  X_train_treat = X,
  treatment_indicator_train = treatment,
  outcome_type = "continuous",
  number_of_trees_control = 50,
  number_of_trees_treat = 50,
  prior_type_control = "horseshoe",
  prior_type_treat = "horseshoe",
  local_hp_control = 0.1 / sqrt(50),
  global_hp_control = 0.1 / sqrt(50),
  local_hp_treat = 0.1 / sqrt(50),
  global_hp_treat = 0.1 / sqrt(50),
  N_post = 1000, N_burn = 500,
  store_posterior_sample = TRUE,
  verbose = TRUE
)

cat("\n")
print(fit_csf)
cat("\n")

tau_hat_csf <- fit_csf$train_predictions_treat
rmse_csf <- sqrt(mean((tau_hat_csf - tau_x)^2))
cat(sprintf("CATE RMSE: %.3f\n", rmse_csf))
cat(sprintf("ATE: %.3f (true: %.3f)\n", mean(tau_hat_csf), mean(tau_x)))

# ── 4. Multi-chain causal ──────────────────────────────────────────────────

cat("\n=== Multi-chain CausalHorseForest (2 chains) ===\n")

set.seed(1)
fit_mc <- CausalHorseForest(
  y = y,
  X_train_control = X,
  X_train_treat = X,
  treatment_indicator_train = treatment,
  outcome_type = "continuous",
  number_of_trees = 50,
  N_post = 500, N_burn = 250,
  store_posterior_sample = TRUE,
  n_chains = 2,
  verbose = FALSE
)
cat(sprintf("Sigma samples: %d (expected %d)\n",
            length(fit_mc$sigma), 500 * 2))
cat(sprintf("Chains: %d\n", fit_mc$mcmc$n_chains))

# ── 5. Causal plots ────────────────────────────────────────────────────────

if (requireNamespace("ggplot2", quietly = TRUE)) {
  cat("\n=== Causal plot tests ===\n")

  p1 <- plot(fit_chf, type = "trace")
  cat("Trace plot OK.\n")

  p2 <- plot(fit_chf, type = "ate")
  cat("ATE plot OK.\n")

  p3 <- plot(fit_chf, type = "cate")
  cat("CATE plot OK.\n")
}

cat("\nAll causal tests passed.\n")
