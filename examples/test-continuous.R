# ── Manual test: continuous outcome models ────────────────────────────────────
#
# Tests HorseTrees and ShrinkageTrees on a continuous outcome with
# known signal, verifying predictions, posterior samples, and multi-chain.
#
# Run from the package root:
#   Rscript tests/manual/test-continuous.R

library(ShrinkageTrees)

# ── 1. Data generation ──────────────────────────────────────────────────────

set.seed(42)
n <- 200; p <- 5
X      <- matrix(rnorm(n * p), ncol = p)
X_test <- matrix(rnorm(50 * p), ncol = p)

# True signal: f(x) = 2*X1 + X2 - 0.5*X3
f_true      <- 2 * X[, 1] + X[, 2] - 0.5 * X[, 3]
f_true_test <- 2 * X_test[, 1] + X_test[, 2] - 0.5 * X_test[, 3]
y <- f_true + rnorm(n, sd = 1)

# ── 2. HorseTrees ───────────────────────────────────────────────────────────

cat("=== HorseTrees (continuous) ===\n")

set.seed(1)
fit_ht <- HorseTrees(
  y = y, X_train = X, X_test = X_test,
  outcome_type = "continuous",
  number_of_trees = 50,
  N_post = 1000, N_burn = 500,
  store_posterior_sample = TRUE,
  verbose = TRUE
)

cat("\n")
print(fit_ht)
cat("\n")
smry <- summary(fit_ht)

# RMSE on train
rmse_train <- sqrt(mean((fit_ht$train_predictions - f_true)^2))
rmse_test  <- sqrt(mean((fit_ht$test_predictions - f_true_test)^2))
cat(sprintf("Train RMSE: %.3f\n", rmse_train))
cat(sprintf("Test  RMSE: %.3f\n", rmse_test))

# Predict on new data
pred <- predict(fit_ht, newdata = X_test)
cat(sprintf("Predict RMSE: %.3f\n", sqrt(mean((pred$mean - f_true_test)^2))))

# ── 3. ShrinkageTrees with different priors ─────────────────────────────────

priors <- c("horseshoe", "standard", "half-cauchy")

cat("\n=== ShrinkageTrees — prior comparison ===\n")
for (pr in priors) {
  set.seed(1)
  fit <- ShrinkageTrees(
    y = y, X_train = X, X_test = X_test,
    outcome_type = "continuous",
    number_of_trees = 50,
    prior_type = pr,
    local_hp  = 0.1 / sqrt(50),
    global_hp = 0.1 / sqrt(50),
    N_post = 500, N_burn = 250,
    verbose = FALSE
  )
  rmse <- sqrt(mean((fit$train_predictions - f_true)^2))
  cat(sprintf("  %-15s  Train RMSE = %.3f  sigma_hat = %.3f\n",
              pr, rmse, mean(fit$sigma)))
}

# ── 4. Multi-chain ──────────────────────────────────────────────────────────

cat("\n=== Multi-chain (2 chains) ===\n")
set.seed(1)
fit_mc <- HorseTrees(
  y = y, X_train = X,
  outcome_type = "continuous",
  number_of_trees = 50,
  N_post = 500, N_burn = 250,
  n_chains = 2,
  verbose = FALSE
)
print(fit_mc)
cat(sprintf("Sigma samples: %d (expected %d)\n",
            length(fit_mc$sigma), 500 * 2))

# ── 5. Plots ────────────────────────────────────────────────────────────────

if (requireNamespace("ggplot2", quietly = TRUE)) {
  cat("\n=== Plot tests ===\n")
  p1 <- plot(fit_mc, type = "trace")
  p2 <- plot(fit_mc, type = "density")
  cat("Trace and density plots created successfully.\n")
}

cat("\nAll continuous tests passed.\n")

