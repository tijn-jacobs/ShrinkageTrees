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
base_seed <- 20250101L
p_fixed   <- 100L
n_trees   <- 200L
n_chains  <- 4L

if (smoke) {
  n_grid <- c(100L, 250L); n_post <- 50L; n_burn <- 50L; n_rep <- 1L
} else {
  n_grid <- c(100L, 250L, 500L, 1000L, 2000L)
  n_post <- 1000L; n_burn <- 1000L; n_rep <- 3L
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
# BART::abart() has no n_chains argument; its multi-chain wrapper mc.abart()
# runs mc.cores independent chains in parallel, mirroring how ShrinkageTrees
# dispatches its own chains via parallel::mclapply.
st_args <- function(d) list(
  y = d$time, status = d$status, X_train = d$X,
  outcome_type = "right-censored", timescale = "time",
  number_of_trees = n_trees, N_post = n_post, N_burn = n_burn,
  n_chains = n_chains, store_posterior_sample = FALSE, verbose = FALSE
)

fits <- list(
  SurvivalBART     = function(d, s) do.call(SurvivalBART, st_args(d)),
  SurvivalDART     = function(d, s) do.call(SurvivalDART, st_args(d)),
  HorseTrees       = function(d, s) do.call(HorseTrees, st_args(d)),
  `BART::mc.abart` = function(d, s) BART::mc.abart(
    x.train = d$X, times = d$time, delta = d$status,
    ntree = n_trees, ndpost = n_post, nskip = n_burn,
    mc.cores = n_chains, seed = s, printevery = 1e6L)
)
if (!has_bart) fits[["BART::mc.abart"]] <- NULL

fn_display <- c(SurvivalBART = "SurvivalBART()", SurvivalDART = "SurvivalDART()",
                HorseTrees = "HorseTrees()", `BART::mc.abart` = "abart()")

time_fit <- function(fn, d, seed) {
  t0  <- Sys.time()
  res <- tryCatch(fn(d, seed), error = function(e) e)
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
agg <- aggregate(seconds ~ fn + n, data = raw[!is.na(raw$seconds), ],
                 FUN = function(z) c(mean = mean(z), sd = sd(z)))
agg <- do.call(data.frame, agg)
names(agg)[3:4] <- c("mean", "sd")
agg$sd[is.na(agg$sd)] <- 0
agg$label <- factor(fn_display[as.character(agg$fn)],
                    levels = unname(fn_display))

if (has_gg) {
  p_scaling <- ggplot2::ggplot(
    agg, ggplot2::aes(x = n, y = mean, colour = label, fill = label)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = pmax(mean - sd, 0), ymax = mean + sd),
      alpha = 0.15, colour = NA) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::geom_point(size = 1.8) +
    ggplot2::labs(x = "Sample size n", y = "Wall-clock seconds",
                  colour = NULL, fill = NULL) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom")
  ggplot2::ggsave("outputs/benchmark_scaling.pdf", p_scaling,
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
