# ── Manual test: coda MCMC diagnostics ────────────────────────────────────────
#
# Tests the as.mcmc.list() S3 method and the automatic convergence diagnostics
# (R-hat, ESS) reported by summary().
#
# Run from the package root:
#   Rscript tests/manual/test-coda-diagnostics.R

library(ShrinkageTrees)

if (!requireNamespace("coda", quietly = TRUE)) {
  stop("This test requires the 'coda' package. ",
       "Install it with: install.packages('coda')")
}

# ── 1. Data generation ──────────────────────────────────────────────────────

set.seed(42)
n <- 200; p <- 5
X <- matrix(rnorm(n * p), ncol = p)
f_true <- 2 * X[, 1] + X[, 2] - 0.5 * X[, 3]
y <- f_true + rnorm(n, sd = 1)

# ── 2. Single-chain fit ──────────────────────────────────────────────────────

cat("=== Single-chain model ===\n")
set.seed(1)
fit_1ch <- HorseTrees(
  y = y, X_train = X,
  outcome_type = "continuous",
  number_of_trees = 50,
  N_post = 500, N_burn = 250,
  verbose = FALSE
)

# 2a. as.mcmc.list -- returns an mcmc.list with 1 chain
mcmc_1ch <- as.mcmc.list(fit_1ch)
cat("Class:", class(mcmc_1ch), "\n")
cat("Number of chains:", length(mcmc_1ch), "\n")
cat("Iterations per chain:", coda::niter(mcmc_1ch), "\n")
cat("Variables:", coda::varnames(mcmc_1ch), "\n")
stopifnot(inherits(mcmc_1ch, "mcmc.list"))
stopifnot(length(mcmc_1ch) == 1)
stopifnot(coda::niter(mcmc_1ch) == 500)

# 2b. ESS (should be available for single chain)
ess_1ch <- coda::effectiveSize(mcmc_1ch)
cat("Effective sample size:", round(ess_1ch), "\n")
stopifnot(ess_1ch > 0 && ess_1ch <= 500)

# 2c. Geweke diagnostic (single chain)
gew <- coda::geweke.diag(mcmc_1ch[[1]])
cat("Geweke z-score:", round(gew$z, 3), "\n")

# 2d. summary() should show ESS but not R-hat (single chain)
smry_1ch <- summary(fit_1ch)
cat("\nSummary diagnostics present:", !is.null(smry_1ch$diagnostics), "\n")
cat("ESS present:", !is.null(smry_1ch$diagnostics$ess), "\n")
cat("R-hat present:", !is.null(smry_1ch$diagnostics$rhat), "\n")
stopifnot(!is.null(smry_1ch$diagnostics$ess))
stopifnot(is.null(smry_1ch$diagnostics$rhat))  # no R-hat for 1 chain

cat("\n--- Single-chain summary ---\n")
print(smry_1ch)

# ── 3. Multi-chain fit ───────────────────────────────────────────────────────

cat("\n=== Multi-chain model (2 chains) ===\n")
set.seed(1)
fit_2ch <- HorseTrees(
  y = y, X_train = X,
  outcome_type = "continuous",
  number_of_trees = 50,
  N_post = 500, N_burn = 250,
  n_chains = 2,
  verbose = FALSE
)

# 3a. as.mcmc.list -- returns mcmc.list with 2 chains
mcmc_2ch <- as.mcmc.list(fit_2ch)
cat("Number of chains:", length(mcmc_2ch), "\n")
cat("Iterations per chain:", coda::niter(mcmc_2ch), "\n")
stopifnot(length(mcmc_2ch) == 2)
stopifnot(coda::niter(mcmc_2ch) == 500)

# 3b. ESS (pooled across chains)
ess_2ch <- coda::effectiveSize(mcmc_2ch)
cat("Effective sample size:", round(ess_2ch), "\n")

# 3c. Gelman-Rubin R-hat
gr <- coda::gelman.diag(mcmc_2ch, autoburnin = FALSE, multivariate = FALSE)
cat("Gelman-Rubin R-hat:\n")
print(gr)
# R-hat should be close to 1 for a well-converged chain
stopifnot(gr$psrf[, "Point est."] < 1.1)

# 3d. summary() should show both ESS and R-hat
smry_2ch <- summary(fit_2ch)
cat("\nR-hat present:", !is.null(smry_2ch$diagnostics$rhat), "\n")
stopifnot(!is.null(smry_2ch$diagnostics$rhat))

cat("\n--- Multi-chain summary ---\n")
print(smry_2ch)

# ── 4. Coda plotting (autocorrelation, crosscorrelation) ─────────────────────

cat("\n=== Coda plots ===\n")

# These are base R plots from coda
cat("Creating autocorrelation plot...\n")
coda::autocorr.plot(mcmc_2ch, auto.layout = TRUE)

cat("Creating traceplot (coda version)...\n")
coda::traceplot(mcmc_2ch)

cat("Creating density plot (coda version)...\n")
coda::densplot(mcmc_2ch)

cat("Creating Gelman plot...\n")
coda::gelman.plot(mcmc_2ch)

# ── 5. Additional coda diagnostics ───────────────────────────────────────────

cat("\n=== Additional diagnostics ===\n")

# Geweke diagnostic per chain
cat("\nGeweke diagnostic (per chain):\n")
for (ch in seq_along(mcmc_2ch)) {
  gew <- coda::geweke.diag(mcmc_2ch[[ch]])
  cat(sprintf("  Chain %d: z = %.3f (p = %.3f)\n",
              ch, gew$z, 2 * pnorm(-abs(gew$z))))
}

# Heidelberger-Welch stationarity test
cat("\nHeidelberger-Welch test:\n")
hw <- coda::heidel.diag(mcmc_2ch)
print(hw)

# Raftery-Lewis diagnostic
cat("\nRaftery-Lewis diagnostic:\n")
rl <- coda::raftery.diag(mcmc_2ch)
print(rl)

# Summary of the mcmc.list itself
cat("\nCoda summary:\n")
print(summary(mcmc_2ch))

# ── 6. Survival model with coda ──────────────────────────────────────────────

cat("\n=== Survival model with coda ===\n")
set.seed(42)
time   <- rexp(n, rate = exp(0.5 * X[, 1]))
status <- rbinom(n, 1, 0.7)

fit_surv <- SurvivalBART(
  time = time, status = status, X_train = X,
  number_of_trees = 10,
  N_post = 300, N_burn = 150,
  n_chains = 2,
  verbose = FALSE
)

mcmc_surv <- as.mcmc.list(fit_surv)
cat("Survival model chains:", length(mcmc_surv), "\n")
cat("Survival model ESS:", round(coda::effectiveSize(mcmc_surv)), "\n")

gr_surv <- coda::gelman.diag(mcmc_surv, autoburnin = FALSE, multivariate = FALSE)
cat("Survival R-hat:", round(gr_surv$psrf[, "Point est."], 3), "\n")

# summary() should include diagnostics for survival model too
smry_surv <- summary(fit_surv)
cat("\n--- Survival model summary ---\n")
print(smry_surv)

# ── 7. Without coda installed (simulated) ────────────────────────────────────

cat("\n=== Verify graceful behaviour ===\n")
# When coda IS installed, diagnostics should be non-NULL
stopifnot(!is.null(smry_2ch$diagnostics))
stopifnot(!is.null(smry_surv$diagnostics))

cat("All assertions passed.\n")

cat("\nAll coda diagnostics tests passed.\n")
