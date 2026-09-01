################################################################################
# Timing benchmark for the ShrinkageTrees R Journal paper.
#
# Reproduces Figure "scaling" (Section 6.5): wall-clock seconds against sample
# size n for SurvivalBART(), SurvivalDART(), HorseTrees() and BART::mc.abart(),
# with p and the MCMC settings held fixed.
#
# Usage
#   Rscript benchmark_timings.R              # full run, settings as in the paper
#   Rscript benchmark_timings.R --smoke      # tiny grid, one replicate
#
# Output (all in outputs/ next to this script)
#   benchmark_timings.rds   raw long-format timings, one row per fit
#   benchmark_timings.csv   the same, as text
#   benchmark_scaling.pdf   the figure used in the manuscript
#
# Timings are hardware dependent. The manuscript reports a 2024 MacBook Air
# (Apple M3, 16 GB) running R 4.5.2, ShrinkageTrees 2.1.0 and BART 2.9.10.
# Absolute seconds will differ elsewhere; the relative ordering of the four
# functions is the claim the figure supports. sessionInfo() is written to the
# .rds so any rerun can be compared against the run behind the figure.
################################################################################

library(ShrinkageTrees)

stopifnot(utils::packageVersion("ShrinkageTrees") >= "2.1.0")
has_bart <- requireNamespace("BART", quietly = TRUE)
has_gg   <- requireNamespace("ggplot2", quietly = TRUE)
if (!has_bart) message("BART not installed: the abart() comparison is skipped.")

args  <- commandArgs(trailingOnly = TRUE)
smoke <- "--smoke" %in% args

# ── Settings (these are the values quoted in Section 6.5) ───────────────────
base_seed <- 2026L
p_fixed   <- 100L
n_trees   <- 200L
n_chains  <- 4L

if (smoke) {
  n_grid <- c(100L, 250L, 500L, 1000L, 2000L)
  n_post <- 100L; n_burn <- 0L; n_rep <- 1L
} else {
  n_grid <- c(100L, 250L, 500L, 1000L, 2000L)
  n_post <- 100L; n_burn <- 0L; n_rep <- 100L
}

# Seed spacing is far wider than any n in the grid, so (replicate, n) pairs
# cannot collide however the grid is extended.
seed_block <- 100000L
stopifnot(max(n_grid) < seed_block)

script_path <- (function() {
  a <- commandArgs(trailingOnly = FALSE)
  m <- grep("--file=", a, fixed = TRUE)
  if (length(m)) sub("--file=", "", a[m], fixed = TRUE) else NULL
})()
if (!is.null(script_path) && nzchar(script_path)) {
  setwd(normalizePath(dirname(script_path)))
}
dir.create("outputs", showWarnings = FALSE)

# fn_display, benchmark_summary, benchmark_figure
source("benchmark_figure.R")

# ── Data-generating process ─────────────────────────────────────────────────
# AFT log-normal, five active covariates out of p, right-censoring at ~35%.
# log T = 2 + x1 - 0.8 x2 + 0.6 x3 - 0.4 x5 + N(0, 0.5^2)
censor_target <- 0.35

calibrate_censoring_rate <- function(target = censor_target, n_ref = 100000L,
                                     seed = 99L) {
  set.seed(seed)
  X     <- matrix(rnorm(n_ref * 5L), n_ref, 5L)
  log_T <- 2 + X[, 1] - 0.8 * X[, 2] + 0.6 * X[, 3] - 0.4 * X[, 5] +
    rnorm(n_ref, sd = 0.5)
  T_ref <- exp(log_T)
  # P(C < T) = E[1 - exp(-rate * T)]; solve for rate.
  uniroot(function(r) mean(1 - exp(-r * T_ref)) - target,
          lower = 1e-8, upper = 100, tol = 1e-6)$root
}

censor_rate <- calibrate_censoring_rate()

sim_data <- function(n, p, seed) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("x", seq_len(p))
  log_T <- 2 + X[, 1] - 0.8 * X[, 2] + 0.6 * X[, 3] - 0.4 * X[, 5] +
    rnorm(n, sd = 0.5)
  T <- exp(log_T)
  C <- rexp(n, rate = censor_rate)
  list(time = pmin(T, C), status = as.integer(T <= C), X = X)
}

# ── Fits ────────────────────────────────────────────────────────────────────
# Three different argument conventions, so three argument builders rather than
# one shared list:
#
#   SurvivalBART()/SurvivalDART() take `time`/`status` and fix `outcome_type`
#     internally, so it must not be supplied. `n_chains` is not a formal
#     argument but reaches ShrinkageTrees() through `...`.
#   HorseTrees() takes `y`/`status` and needs `outcome_type` explicitly.
#   BART::mc.abart() takes `times`/`delta`.
wrapper_args <- function(d) list(
  time = d$time, status = d$status, X_train = d$X,
  timescale = "time",
  number_of_trees = n_trees, N_post = n_post, N_burn = n_burn,
  n_chains = n_chains, store_posterior_sample = FALSE, verbose = FALSE
)

horsetrees_args <- function(d) list(
  y = d$time, status = d$status, X_train = d$X,
  outcome_type = "right-censored", timescale = "time",
  number_of_trees = n_trees, N_post = n_post, N_burn = n_burn,
  n_chains = n_chains, store_posterior_sample = FALSE, verbose = FALSE
)

# BART::mc.abart() splits `ndpost` ACROSS cores: with ndpost = 50 and
# mc.cores = 4 it keeps ceiling(50/4) = 13 draws per chain. ShrinkageTrees
# treats N_post as per chain. Multiplying here is what makes the comparison
# apples to apples; without it abart does a quarter of the posterior work and
# looks correspondingly faster. `nskip` is per chain in both.
fits <- list(
  SurvivalBART     = function(d, s) do.call(SurvivalBART, wrapper_args(d)),
  SurvivalDART     = function(d, s) do.call(SurvivalDART, wrapper_args(d)),
  HorseTrees       = function(d, s) do.call(HorseTrees, horsetrees_args(d)),
  `BART::mc.abart` = function(d, s) BART::mc.abart(
    x.train = d$X, times = d$time, delta = d$status,
    ntree = n_trees, ndpost = n_post * n_chains, nskip = n_burn,
    mc.cores = n_chains, seed = s, printevery = 1e6L)
)
if (!has_bart) fits[["BART::mc.abart"]] <- NULL

time_fit <- function(fn, d, seed) {
  # BART prints a data/prior header regardless of printevery; swallow it so
  # the timing log stays readable. The capture is outside the clock.
  t0  <- Sys.time()
  res <- tryCatch(utils::capture.output(val <- fn(d, seed)),
                  error = function(e) e)
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  list(seconds = if (inherits(res, "error")) NA_real_ else secs,
       error   = if (inherits(res, "error")) conditionMessage(res) else NA_character_)
}

# ── Run ─────────────────────────────────────────────────────────────────────
cat("TIMING BENCHMARK\n")
cat("  ShrinkageTrees ", format(utils::packageVersion("ShrinkageTrees")),
    if (smoke) " | SMOKE RUN" else "", "\n", sep = "")
cat("  n in {", paste(n_grid, collapse = ", "), "} | p = ", p_fixed,
    " | m = ", n_trees, " | N_post = ", n_post, " | n_chains = ", n_chains,
    " | ", n_rep, " replicates\n\n", sep = "")

rows <- list()
for (n in n_grid) {
  for (r in seq_len(n_rep)) {
    seed <- base_seed + seed_block * r + n   # data shared by all functions
    d    <- sim_data(n, p_fixed, seed)
    for (nm in names(fits)) {
      tm <- time_fit(fits[[nm]], d, seed)
      cat(sprintf("  n = %5d  rep %d  %-16s %8.2f s%s\n", n, r, nm,
                  tm$seconds, if (is.na(tm$error)) "" else "  FAILED"))
      rows[[length(rows) + 1L]] <- data.frame(
        fn = nm, n = n, p = p_fixed, replicate = r, seed = seed,
        seconds = tm$seconds, error = tm$error, stringsAsFactors = FALSE)
    }
  }
}

raw <- do.call(rbind, rows)
raw$fn <- factor(raw$fn, levels = names(fn_display)[names(fn_display) %in% raw$fn])

attr(raw, "session_info")    <- utils::sessionInfo()
attr(raw, "censoring_rate")  <- censor_rate
saveRDS(raw, "outputs/benchmark_timings.rds")
write.csv(raw, "outputs/benchmark_timings.csv", row.names = FALSE)

# ── Figure ──────────────────────────────────────────────────────────────────
agg <- benchmark_summary(raw)

if (has_gg) {
  ggplot2::ggsave("outputs/benchmark_scaling.pdf", benchmark_figure(agg),
                  width = 6, height = 3.5)
  cat("\nWrote outputs/benchmark_scaling.pdf\n")
} else {
  cat("\nggplot2 not installed: figure skipped, timings still written.\n")
}

print(agg[, c("label", "n", "mean", "sd")])

failed <- raw[!is.na(raw$error), ]
if (nrow(failed)) {
  cat("\nFAILURES:\n"); print(failed[, c("fn", "n", "replicate", "error")])
}
