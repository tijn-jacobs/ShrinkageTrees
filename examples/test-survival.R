# ── Manual test: survival outcome models ──────────────────────────────────────
#
# Tests HorseTrees with right-censored and interval-censored outcomes.
# Also tests survival wrappers (SurvivalBART, SurvivalBCF).
#
# Run from the package root:
#   Rscript tests/manual/test-survival.R

library(ShrinkageTrees)

# ── 1. Right-censored survival ──────────────────────────────────────────────

cat("=== Right-censored survival ===\n")

set.seed(42)
n <- 200; p <- 5
X      <- matrix(rnorm(n * p), ncol = p)
X_test <- matrix(rnorm(50 * p), ncol = p)

# Generate log-normal survival times
linpred <- X[, 1] - 0.5 * X[, 2]
true_time <- exp(linpred + rnorm(n, sd = 0.5))

# Random censoring (~30%)
censor_time <- rexp(n, rate = 0.3)
time_obs <- pmin(true_time, censor_time)
status   <- as.numeric(true_time <= censor_time)
cat(sprintf("Censoring rate: %.1f%%\n", 100 * mean(status == 0)))

set.seed(1)
fit_rc <- HorseTrees(
  y = time_obs, status = status,
  X_train = X, X_test = X_test,
  outcome_type = "right-censored",
  timescale = "time",
  number_of_trees = 50,
  N_post = 500, N_burn = 250,
  store_posterior_sample = TRUE,
  verbose = TRUE
)

cat("\n")
print(fit_rc)
cat("\n")
summary(fit_rc)

# Predictions should be positive
stopifnot(all(fit_rc$train_predictions > 0))
stopifnot(all(fit_rc$test_predictions > 0))

# C-index
if (requireNamespace("survival", quietly = TRUE)) {
  cindex <- survival::concordance(
    survival::Surv(time_obs, status) ~ fit_rc$train_predictions
  )$concordance
  cat(sprintf("C-index: %.3f\n", cindex))
}

# ── 2. Interval-censored survival ──────────────────────────────────────────

cat("\n=== Interval-censored survival ===\n")

set.seed(42)
true_time2 <- exp(X[, 1] + rnorm(n, sd = 0.5))

# Mixed censoring: 40% exact, 30% interval, 30% right-censored
obs_type <- sample(c("exact", "interval", "right"),
                   n, replace = TRUE, prob = c(0.4, 0.3, 0.3))
left_time  <- numeric(n)
right_time <- numeric(n)

for (i in seq_len(n)) {
  if (obs_type[i] == "exact") {
    left_time[i]  <- true_time2[i]
    right_time[i] <- true_time2[i]
  } else if (obs_type[i] == "interval") {
    w <- runif(1, 0.1, 0.5) * true_time2[i]
    left_time[i]  <- true_time2[i] - w / 2
    right_time[i] <- true_time2[i] + w / 2
  } else {
    left_time[i]  <- true_time2[i] * runif(1, 0.3, 0.9)
    right_time[i] <- Inf
  }
}

set.seed(1)
fit_ic <- HorseTrees(
  left_time = left_time, right_time = right_time,
  X_train = X, X_test = X_test,
  outcome_type = "interval-censored",
  timescale = "time",
  number_of_trees = 50,
  N_post = 500, N_burn = 250,
  store_posterior_sample = TRUE,
  verbose = TRUE
)

cat("\n")
print(fit_ic)

stopifnot(all(fit_ic$train_predictions > 0))
stopifnot(all(fit_ic$test_predictions > 0))
cat("Interval-censored model OK.\n")

# ── 3. Survival wrappers ───────────────────────────────────────────────────

cat("\n=== SurvivalBART wrapper ===\n")
set.seed(1)
fit_sbart <- SurvivalBART(
  time = time_obs, status = status,
  X_train = X, X_test = X_test,
  number_of_trees = 50,
  N_post = 500, N_burn = 250,
  verbose = FALSE
)
print(fit_sbart)
stopifnot(all(fit_sbart$train_predictions > 0))
cat("SurvivalBART OK.\n")

# ── 4. Survival plots ─────────────────────────────────────────────────────

if (requireNamespace("ggplot2", quietly = TRUE)) {
  cat("\n=== Survival plot tests ===\n")

  # Population-averaged survival curve
  p1 <- plot(fit_rc, type = "survival")
  cat("Population-averaged survival curve OK.\n")
  p1

  # Individual curves
  p2 <- plot(fit_rc, type = "survival", obs = c(1, 5, 10))
  cat("Individual survival curves OK.\n")
  p2

  # KM overlay
  if (requireNamespace("survival", quietly = TRUE)) {
    p3 <- plot(fit_rc, type = "survival", km = TRUE)
    cat("KM overlay OK.\n")
  }
  p3
}

cat("\nAll survival tests passed.\n")

