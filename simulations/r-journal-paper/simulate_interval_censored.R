################################################################################
# Interval-censored simulation for the ShrinkageTrees R Journal paper.
#
# Reproduces Table "sim-ic-table" (Section 4.4) in full: three priors
# (SurvivalBART, SurvivalDART, HorseTrees) x three dimensions p in {50, 500,
# 5000} x 1000 replicates, evaluated on train and test samples.
#
# Usage
#   Rscript simulate_interval_censored.R              # full run
#   Rscript simulate_interval_censored.R --p 500      # one dimension
#   Rscript simulate_interval_censored.R --cores 8
#   Rscript simulate_interval_censored.R --out /scratch/me/ic
#   Rscript simulate_interval_censored.R 192            # bare core count
#
# On a cluster the core count is taken from the scheduler
# (SLURM_CPUS_PER_TASK and friends) and output goes to $TMPDIR unless --out is
# given. $TMPDIR is node-local and wiped at the end of the job, so copy the
# .rds files somewhere permanent before the job exits.
#
# Output (in $TMPDIR, --out, or outputs/ next to this script)
#   simulate_interval_censored_output.rds  all dimensions combined
#
# Reproducibility
#   Replicate i at dimension p uses seed base_seed + seed_block * b + i, where
#   b is the position of p in known_p. The seed is set inside the worker, so
#   results do not depend on the number of cores, on the scheduling order, or
#   on whether the dimensions are run together or separately: --p 500 alone
#   reproduces exactly the p = 500 rows of a full run.
#
# Requires ShrinkageTrees >= 2.1.0. Earlier versions contain two defects in the
# horseshoe global update and a different default for k; see NEWS.md.
################################################################################

library(ShrinkageTrees)
library(doParallel)
library(foreach)

stopifnot(utils::packageVersion("ShrinkageTrees") >= "2.1.0")

# ── Command line ────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
opt  <- function(flag, default) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) default else args[i + 1L]
}
p_arg   <- opt("--p", NA_character_)

# detectCores() reports the machine, not the allocation, so on a cluster it
# oversubscribes badly. Prefer what the scheduler granted.
resolve_cores <- function() {
  from_flag <- opt("--cores", NA_character_)
  if (!is.na(from_flag)) return(as.integer(from_flag))
  # Bare positional count, i.e. `Rscript simulate_interval_censored.R 192`.
  if (length(args) >= 1L && !grepl("^--", args[1L])) {
    v <- suppressWarnings(as.integer(args[1L]))
    if (!is.na(v)) return(v)
  }
  # Deliberately NOT OMP_NUM_THREADS: job scripts set it to 1 to stop BLAS
  # threading inside each worker, which is unrelated to how many workers to run.
  for (v in c("SLURM_CPUS_PER_TASK", "SLURM_CPUS_ON_NODE", "NSLOTS", "PBS_NP")) {
    val <- Sys.getenv(v)
    if (nzchar(val) && !is.na(suppressWarnings(as.integer(val))))
      return(as.integer(val))
  }
  max(1L, parallel::detectCores() - 1L)   # local: leave one core free
}
n_cores <- resolve_cores()

# ── Settings ────────────────────────────────────────────────────────────────
base_seed  <- 2026L
seed_block <- 1000000L          # far wider than any plausible n_rep
known_p    <- c(50L, 500L, 5000L)

p_values <- if (is.na(p_arg)) known_p else as.integer(strsplit(p_arg, ",")[[1]])
if (any(is.na(match(p_values, known_p))))
  stop("Every p must appear in known_p, which fixes its seed block. ",
       "Add the new dimension to known_p rather than relying on its position ",
       "in p_values, so that existing dimensions keep their seeds.")

n_rep <- 1000L; n_obs <- 200L; n_test <- 1000L
n_post <- 5000L; n_burn <- 5000L; n_trees <- 200L
n_chains <- 1L
k_horse  <- 1.0        # 2.1.0 default; the old 0.1 was calibrated pre-fix

stopifnot(n_rep < seed_block)   # guarantees the seed blocks cannot overlap

# Run from anywhere: work relative to this script.
script_path <- (function() {
  a <- commandArgs(trailingOnly = FALSE)
  m <- grep("--file=", a, fixed = TRUE)
  if (length(m)) sub("--file=", "", a[m], fixed = TRUE) else NULL
})()
if (!is.null(script_path) && nzchar(script_path)) {
  setwd(normalizePath(dirname(script_path)))
}

# --out wins; then TMPDIR, which cluster jobs set and which is node-local and
# fast; then outputs/ beside the script. Copy results off TMPDIR before the
# job ends, since it is wiped.
out_dir <- opt("--out", NA_character_)
if (is.na(out_dir)) {
  tmp <- Sys.getenv("TMPDIR")
  out_dir <- if (nzchar(tmp)) tmp else "outputs"
}
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ── Data-generating process ─────────────────────────────────────────────────
# Friedman + sparse linear (spike-and-slab) on X ~ U[0,1]^p. beta is drawn
# fresh in each call, so every replicate sees a new active set. f is rescaled
# to unit variance via a reference draw, and sigma is set so that
# var(Y) / sigma^2 = var_ratio, i.e. sigma^2 = 1 / (var_ratio - 1).
data_gen <- function(n, n_test, p,
                     var_ratio = 10 / 9, s_slab = 0.05, n_ref = 5000L) {
  beta         <- rep(0, p)
  active       <- rbinom(p, 1L, s_slab) == 1L
  beta[active] <- rnorm(sum(active))

  f_raw <- function(x) {
    10 * sin(pi * x[, 1] * x[, 2]) +
      20 * (x[, 3] - 0.5)^2 +
      10 * x[, 4] +
       5 * x[, 5] +
      as.vector(x %*% beta)
  }

  x_ref <- matrix(runif(n_ref * p), n_ref, p)
  sd_f  <- sd(f_raw(x_ref))
  rm(x_ref)

  f     <- function(x) f_raw(x) / sd_f
  sigma <- sqrt(1 / (var_ratio - 1))

  x_train <- matrix(runif(n      * p), n,      p)
  x_test  <- matrix(runif(n_test * p), n_test, p)
  colnames(x_train) <- colnames(x_test) <- paste0("x", seq_len(p))

  f_train <- f(x_train)
  f_test  <- f(x_test)
  log_t   <- f_train + rnorm(n, 0, sigma)

  # Interval censoring on the log scale: three inspection times per subject.
  log_t_floor <- min(log_t) - 1
  log_t_max   <- quantile(log_t, 0.95)
  v <- t(apply(matrix(runif(n * 3, log_t_floor, log_t_max), n, 3), 1, sort))

  log_left  <- numeric(n)
  log_right <- numeric(n)
  for (i in seq_len(n)) {
    li <- log_t[i]
    if (li > v[i, 3]) {
      log_left[i]  <- v[i, 3]
      log_right[i] <- Inf
    } else if (li <= v[i, 1]) {
      log_left[i]  <- log_t_floor
      log_right[i] <- v[i, 1]
    } else {
      j <- which(v[i, ] >= li)[1]
      log_left[i]  <- v[i, j - 1]
      log_right[i] <- v[i, j]
    }
  }

  list(x_train = x_train, x_test = x_test,
       f_train = f_train, f_test = f_test,
       log_left = log_left, log_right = log_right,
       sigma = sigma, n_active = sum(active),
       n_rcens = sum(!is.finite(log_right)))
}

# ── Metrics ─────────────────────────────────────────────────────────────────
compute_metrics <- function(f_hat, f_post, f_true) {
  ci_lower <- apply(f_post, 2, quantile, 0.025)
  ci_upper <- apply(f_post, 2, quantile, 0.975)
  c(RMSE     = sqrt(mean((f_hat - f_true)^2)),
    Coverage = mean(f_true >= ci_lower & f_true <= ci_upper),
    Length   = mean(ci_upper - ci_lower))
}

# ── Fitting ─────────────────────────────────────────────────────────────────
shared_args <- function(data, n_post, n_burn, n_chains, n_trees) list(
  left_time  = data$log_left,
  right_time = data$log_right,
  timescale  = "log",
  X_train    = data$x_train,
  X_test     = data$x_test,
  number_of_trees = n_trees,
  N_post     = n_post,
  N_burn     = n_burn,
  n_chains   = n_chains,
  store_posterior_sample = TRUE,
  verbose    = FALSE
)

fit_outputs <- function(fit) list(
  train = list(f_hat  = as.numeric(fit$train_predictions),
               f_post = fit$train_predictions_sample),
  test  = list(f_hat  = as.numeric(fit$test_predictions),
               f_post = fit$test_predictions_sample)
)

methods <- c("BART", "DART", "Horseshoe Forest")

fit_method <- function(method, data, n_post, n_burn, n_chains, n_trees, k) {
  a <- shared_args(data, n_post, n_burn, n_chains, n_trees)
  fit_outputs(switch(
    method,
    BART               = do.call(SurvivalBART, a),
    DART               = do.call(SurvivalDART, a),
    `Horseshoe Forest` = do.call(HorseTrees,
                                 c(a, list(outcome_type = "interval-censored",
                                           k = k)))
  ))
}

# ── One replicate ───────────────────────────────────────────────────────────
run_one_sim <- function(seed, n, n_test, p, n_post, n_burn, n_chains,
                        n_trees, k) {
  set.seed(seed)
  dat <- data_gen(n = n, n_test = n_test, p = p)

  out <- do.call(rbind, lapply(methods, function(method) {
    fit <- fit_method(method, dat, n_post, n_burn, n_chains, n_trees, k)
    do.call(rbind, lapply(c("train", "test"), function(split) {
      truth <- if (split == "train") dat$f_train else dat$f_test
      mv <- compute_metrics(fit[[split]]$f_hat, fit[[split]]$f_post, truth)
      data.frame(Method = method, Split = split,
                 RMSE = unname(mv["RMSE"]),
                 Coverage = unname(mv["Coverage"]),
                 Length = unname(mv["Length"]),
                 stringsAsFactors = FALSE)
    }))
  }))

  out$n_active <- dat$n_active
  out$n_rcens  <- dat$n_rcens
  out$seed     <- seed
  out
}

# ── Run ─────────────────────────────────────────────────────────────────────
registerDoParallel(cores = n_cores)
cat("SIMULATION: interval-censored (Friedman + spike-and-slab)\n")
cat("  ShrinkageTrees ", format(utils::packageVersion("ShrinkageTrees")),
    " | cores ", n_cores, "\n", sep = "")
cat("  output -> ", normalizePath(out_dir), "\n", sep = "")
cat("  p in {", paste(p_values, collapse = ", "), "} | ", n_rep,
    " replicates | n = ", n_obs, ", n_test = ", n_test,
    " | m = ", n_trees, ", N_post = ", n_post, ", k = ", k_horse,
    "\n\n", sep = "")

all_res <- list()
for (p_val in p_values) {
  # The seed block is keyed to the value of p, not to its position in
  # p_values, so a single-dimension rerun reproduces that dimension exactly.
  offset <- base_seed + seed_block * match(p_val, known_p)
  cat("p = ", p_val, " (seeds ", offset + 1L, " .. ", offset + n_rep, ")\n",
      sep = "")
  t0 <- Sys.time()

  res <- foreach(i = seq_len(n_rep), .combine = "rbind",
                 .packages = "ShrinkageTrees") %dopar% {
    r <- run_one_sim(seed = offset + i, n = n_obs, n_test = n_test, p = p_val,
                     n_post = n_post, n_burn = n_burn, n_chains = n_chains,
                     n_trees = n_trees, k = k_horse)
    r$Iter <- i
    r$p    <- p_val
    r
  }

  cat("  done in ", format(round(difftime(Sys.time(), t0), 1)), "\n", sep = "")
  all_res[[as.character(p_val)]] <- res
}

combined <- do.call(rbind, all_res)
combined$Method <- factor(combined$Method, levels = methods)
combined$Split  <- factor(combined$Split,  levels = c("train", "test"))
attr(combined, "session_info") <- utils::sessionInfo()
saveRDS(combined, file.path(out_dir, "simulate_interval_censored_output.rds"))

cat("\nWrote ", file.path(out_dir, "simulate_interval_censored_output.rds"),
    " (", nrow(combined), " rows)\n", sep = "")
print(aggregate(cbind(RMSE, Coverage, Length) ~ p + Method + Split,
                data = combined, FUN = mean))
