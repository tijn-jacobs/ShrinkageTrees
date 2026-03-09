# ── Manual test: binary outcome ───────────────────────────────────────────────
#
# Tests HorseTrees with a binary outcome (probit BART).
#
# Run from the package root:
#   Rscript tests/manual/test-binary.R

library(ShrinkageTrees)

# ── 1. Data generation ──────────────────────────────────────────────────────

set.seed(42)
n <- 200; p <- 5
X      <- matrix(rnorm(n * p), ncol = p)
X_test <- matrix(rnorm(50 * p), ncol = p)

# True latent: X1 + 0.5*X2
latent <- X[, 1] + 0.5 * X[, 2]
prob_true <- pnorm(latent)
y <- rbinom(n, 1, prob_true)

cat(sprintf("Prevalence: %.1f%%\n", 100 * mean(y)))

# ── 2. Fit ──────────────────────────────────────────────────────────────────

cat("\n=== HorseTrees (binary) ===\n")

set.seed(1)
fit <- HorseTrees(
  y = y, X_train = X, X_test = X_test,
  outcome_type = "binary",
  number_of_trees = 50,
  N_post = 1000, N_burn = 500,
  store_posterior_sample = TRUE,
  verbose = TRUE
)

cat("\n")
print(fit)
cat("\n")
summary(fit)

# Probabilities should be in [0, 1]
stopifnot(all(fit$train_probabilities >= 0 & fit$train_probabilities <= 1))
stopifnot(all(fit$test_probabilities  >= 0 & fit$test_probabilities  <= 1))

# Brier score
brier <- mean((fit$train_probabilities - y)^2)
cat(sprintf("Brier score: %.4f\n", brier))

# AUC (manual, no extra dependency)
pos <- fit$train_probabilities[y == 1]
neg <- fit$train_probabilities[y == 0]
auc <- mean(outer(pos, neg, ">")) + 0.5 * mean(outer(pos, neg, "=="))
cat(sprintf("AUC: %.3f\n", auc))

# Predict on new data
pred <- predict(fit, newdata = X_test)
stopifnot(all(pred$probabilities >= 0 & pred$probabilities <= 1))
cat(sprintf("Test predictions: %d observations\n", length(pred$mean)))

cat("\nAll binary tests passed.\n")

