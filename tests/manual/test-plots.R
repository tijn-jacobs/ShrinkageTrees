# ── Manual test: plot methods ─────────────────────────────────────────────────
#
# Visual verification of all plot types for ShrinkageTrees and
# CausalShrinkageForest. Saves plots to a temporary PDF so they can be
# reviewed.
#
# Run from the package root:
#   Rscript tests/manual/test-plots.R

library(ShrinkageTrees)

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("ggplot2 is required for plot tests.")
}

outfile <- tempfile("plot-tests-", fileext = ".pdf")
cat("Saving plots to:", outfile, "\n\n")
pdf(outfile, width = 8, height = 5)

# ── 1. ShrinkageTrees plots ────────────────────────────────────────────────

set.seed(1)
n <- 100; p <- 5
X <- matrix(rnorm(n * p), ncol = p)
y <- X[, 1] + rnorm(n)

fit_ht <- HorseTrees(
  y = y, X_train = X,
  outcome_type = "continuous",
  number_of_trees = 20,
  N_post = 200, N_burn = 100,
  store_posterior_sample = TRUE,
  n_chains = 2,
  verbose = FALSE
)

cat("1. Sigma traceplot (2 chains)...\n")
print(plot(fit_ht, type = "trace"))

cat("2. Sigma density (2 chains)...\n")
print(plot(fit_ht, type = "density"))

# Dirichlet model for VI plot
fit_dart <- ShrinkageTrees(
  y = y, X_train = X,
  outcome_type = "continuous",
  number_of_trees = 20,
  prior_type = "dirichlet",
  local_hp = 1 / sqrt(20),
  N_post = 200, N_burn = 100,
  verbose = FALSE
)

cat("3. Variable importance (Dirichlet)...\n")
print(plot(fit_dart, type = "vi", n_vi = 5))

# Survival plot
time <- rexp(n, rate = exp(0.5 * X[, 1]))
status <- rbinom(n, 1, 0.7)

fit_surv <- SurvivalBART(
  time = time, status = status, X_train = X,
  number_of_trees = 20, N_post = 200, N_burn = 100,
  store_posterior_sample = TRUE,
  verbose = FALSE
)

cat("4. Population-averaged survival curve...\n")
print(plot(fit_surv, type = "survival"))

cat("5. Individual survival curves...\n")
print(plot(fit_surv, type = "survival", obs = c(1, 5, 10)))

if (requireNamespace("survival", quietly = TRUE)) {
  cat("6. Survival curve with KM overlay...\n")
  print(plot(fit_surv, type = "survival", km = TRUE))
}

# ── 2. CausalShrinkageForest plots ─────────────────────────────────────────

treatment <- rbinom(n, 1, 0.5)
y_causal  <- X[, 1] + (treatment - 0.5) * (1 + X[, 2]) + rnorm(n, sd = 0.5)

fit_causal <- CausalHorseForest(
  y = y_causal,
  X_train_control = X, X_train_treat = X,
  treatment_indicator_train = treatment,
  outcome_type = "continuous",
  number_of_trees = 20,
  N_post = 200, N_burn = 100,
  store_posterior_sample = TRUE,
  n_chains = 2,
  verbose = FALSE
)

cat("7. Causal sigma traceplot (2 chains)...\n")
print(plot(fit_causal, type = "trace"))

cat("8. Causal sigma density (2 chains)...\n")
print(plot(fit_causal, type = "density"))

cat("9. ATE posterior...\n")
print(plot(fit_causal, type = "ate"))

cat("10. CATE caterpillar plot...\n")
print(plot(fit_causal, type = "cate"))

dev.off()
cat("\nAll plots saved to:", outfile, "\n")
cat("Open the PDF to visually inspect.\n")
