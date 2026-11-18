#' Resolve the number of worker cores for multi-chain MCMC
#'
#' Returns a safe core count that never exceeds CRAN's two-core policy
#' when running under `R CMD check`, honours the `OMP_THREAD_LIMIT`
#' environment variable, and falls back to sequential execution on Windows
#' (where `mclapply` cannot fork).
#'
#' @param n_chains Requested number of chains.
#' @return Integer number of cores to use.
#' @keywords internal
#' @noRd
.resolve_cores <- function(n_chains) {
  if (.Platform$OS.type == "windows") return(1L)
  cap <- if (nzchar(Sys.getenv("_R_CHECK_LIMIT_CORES_"))) 2L else Inf
  tl <- suppressWarnings(as.integer(Sys.getenv("OMP_THREAD_LIMIT", "")))
  if (!is.na(tl) && tl > 0) cap <- min(cap, tl)
  detected <- tryCatch(parallel::detectCores(logical = FALSE),
                       error = function(e) 1L)
  if (is.na(detected) || detected < 1L) detected <- 1L
  as.integer(max(1L, min(n_chains, detected, cap)))
}


#' Derive status, censoring indicator and starting values from interval bounds
#'
#' Shared by the fitting functions and the survival wrappers, so that the
#' encoding of interval-censored data is defined in one place.
#'
#' @param left_time,right_time Numeric vectors of interval bounds. An exact
#'   event has `left_time == right_time`; a right-censored observation has
#'   `right_time = Inf`.
#'
#' @return A list with elements:
#'   \item{status}{1 for an exact event, 0 otherwise.}
#'   \item{ic_indicator}{1 for a finite interval with `left < right`, 0
#'     otherwise.}
#'   \item{y}{Starting value for the latent event time: the event time, the
#'     interval midpoint, or the censoring time.}
#' @noRd
.interval_censoring <- function(left_time, right_time) {
  status <- as.integer(left_time == right_time)
  ic_indicator <- as.integer(left_time < right_time & is.finite(right_time))
  y <- ifelse(status == 1, left_time,
              ifelse(ic_indicator == 1, (left_time + right_time) / 2,
                     left_time))
  list(status = status, ic_indicator = ic_indicator, y = y)
}


#' Compute mean estimate for censored data
#'
#' Estimates the mean and standard deviation for right-censored or
#' interval-censored survival data. Uses the `survival` package if available,
#' otherwise falls back to the naive mean among observed events.
#'
#' @param y Numeric vector of (log-transformed) survival times.
#' @param status Numeric vector; event indicator (1 = event, 0 = censored).
#' @param left_time Optional numeric vector of left boundaries (for interval censoring).
#' @param right_time Optional numeric vector of right boundaries (for interval censoring).
#'
#' @return A list with elements:
#'   \item{mu}{Estimated mean of survival times.}
#'   \item{sd}{Estimated standard deviation of survival times.}
#'
#' @noRd
censored_info <- function(y, status, left_time = NULL, right_time = NULL) {

  if (requireNamespace("survival", quietly = TRUE)) {

    if (!is.null(left_time) && !is.null(right_time)) {
      # Interval-censored path: use Surv(time, time2, type = "interval2")
      #
      # Under interval2 semantics an open right end must be NA (or Inf); a
      # finite right_time equal to left_time means an EXACT event. Callers must
      # therefore pass the true bounds here, before any collapse of Inf onto
      # left_time for the C++ layer. Normalising non-finite bounds to NA makes
      # that explicit and is what survreg() expects.
      right_time[!is.finite(right_time)] <- NA_real_
      surv_obj <- survival::Surv(left_time, right_time, type = "interval2")
    } else {
      # Right-censored path
      surv_obj <- survival::Surv(y, status)
    }

    fit <- survival::survreg(surv_obj ~ 1, dist = "gaussian")
    mu <- as.numeric(fit$coefficients[1])
    sd <- as.numeric(fit$scale)

  } else {

    # Fallback: use only uncensored observations
    mu <- mean(y[status == 1], na.rm = TRUE)
    sd <- sd(y[status == 1], na.rm = TRUE)
  }

  list(mu = mu, sd = sd)
}


#' Run several independent chains of a fitting function
#'
#' @param fun The fitting function to call once per chain.
#' @param args Named list of evaluated arguments for `fun`. `n_chains` and
#'   `verbose` are overridden so each call runs a single silent chain.
#' @param n_chains Number of chains.
#' @param verbose Whether to print a one-line progress message.
#'
#' @return A list of fitted objects, one per chain.
#' @noRd
.run_chains <- function(fun, args, n_chains, verbose) {
  args$n_chains <- 1
  args$verbose  <- FALSE

  n_cores <- .resolve_cores(n_chains)
  if (verbose)
    message("Running ", n_chains, " chains (", n_cores, " cores) ...")

  parallel::mclapply(
    seq_len(n_chains),
    function(i) do.call(fun, args),
    mc.cores = n_cores
  )
}

# Row-bind (matrices) or concatenate (vectors) one field across chains.
# Both return NULL when no chain has the field.
.pool_matrix <- function(chains, field) {
  mats <- lapply(chains, `[[`, field)
  if (all(vapply(mats, is.null, logical(1)))) return(NULL)
  do.call(rbind, mats)
}
.pool_vector <- function(chains, field) {
  vecs <- lapply(chains, `[[`, field)
  if (all(vapply(vecs, is.null, logical(1)))) return(NULL)
  unlist(vecs)
}

# Internal helper: pool a list of ShrinkageTrees objects (one per chain) into
# a single ShrinkageTrees object with concatenated/averaged posteriors.
.combine_chains <- function(chains) {
  n_chains <- length(chains)
  out      <- chains[[1]]   # use first chain as structural template

  pool_matrix <- function(field) .pool_matrix(chains, field)
  pool_vector <- function(field) .pool_vector(chains, field)

  # ── Posterior sample matrices (rows = MCMC iterations) ───────────────────
  out$train_predictions_sample <- pool_matrix("train_predictions_sample")
  out$test_predictions_sample  <- pool_matrix("test_predictions_sample")

  # ── Posterior means: recompute from pooled samples when available ─────────
  if (!is.null(out$train_predictions_sample)) {
    out$train_predictions <- colMeans(out$train_predictions_sample)
    out$test_predictions  <- colMeans(out$test_predictions_sample)
  } else {
    # store_posterior_sample = FALSE: average the per-chain posterior means
    out$train_predictions <- Reduce("+", lapply(chains, `[[`, "train_predictions")) / n_chains
    out$test_predictions  <- Reduce("+", lapply(chains, `[[`, "test_predictions"))  / n_chains
  }

  # ── Sigma posterior samples ───────────────────────────────────────────────
  out$sigma <- pool_vector("sigma")

  # ── Acceptance ratios (callers always use mean(); we concatenate) ─────────
  out$acceptance_ratio <- pool_vector("acceptance_ratio")

  # ── Binary outcome: pool probability samples ──────────────────────────────
  out$train_probabilities_sample <- pool_matrix("train_probabilities_sample")
  out$test_probabilities_sample  <- pool_matrix("test_probabilities_sample")
  if (!is.null(out$train_probabilities_sample)) {
    out$train_probabilities <- colMeans(out$train_probabilities_sample)
    out$test_probabilities  <- colMeans(out$test_probabilities_sample)
  } else if (!is.null(chains[[1]]$train_probabilities)) {
    out$train_probabilities <- Reduce("+", lapply(chains, `[[`, "train_probabilities")) / n_chains
    out$test_probabilities  <- Reduce("+", lapply(chains, `[[`, "test_probabilities"))  / n_chains
  }

  # ── Variable importance: rbind iteration matrices so vi helpers still work ─
  out$split_probs        <- pool_matrix("split_probs")
  out$store_split_counts <- pool_matrix("store_split_counts")

  # ── Chain-level diagnostics ───────────────────────────────────────────────
  out$chains <- list(
    n_chains         = n_chains,
    acceptance_ratios = vapply(chains, function(ch) mean(ch$acceptance_ratio), numeric(1))
  )

  # Record n_chains in mcmc slot (N_post stays per-chain for predict())
  out$mcmc$n_chains <- n_chains

  out
}

# Internal helper: pool a list of CausalShrinkageForest objects (one per chain)
# into a single object with concatenated/averaged posteriors.
.combine_causal_chains <- function(chains) {
  n_chains <- length(chains)
  out      <- chains[[1]]

  pool_matrix <- function(field) .pool_matrix(chains, field)
  pool_vector <- function(field) .pool_vector(chains, field)
  avg <- function(field) Reduce("+", lapply(chains, `[[`, field)) / n_chains

  # ── Posterior sample matrices ─────────────────────────────────────────────
  out$train_predictions_sample_control <- pool_matrix("train_predictions_sample_control")
  out$test_predictions_sample_control  <- pool_matrix("test_predictions_sample_control")
  out$train_predictions_sample_treat   <- pool_matrix("train_predictions_sample_treat")
  out$test_predictions_sample_treat    <- pool_matrix("test_predictions_sample_treat")

  # ── Posterior means ───────────────────────────────────────────────────────
  if (!is.null(out$train_predictions_sample_control)) {
    out$train_predictions_control <- colMeans(out$train_predictions_sample_control)
    out$test_predictions_control  <- colMeans(out$test_predictions_sample_control)
    out$train_predictions_treat   <- colMeans(out$train_predictions_sample_treat)
    out$test_predictions_treat    <- colMeans(out$test_predictions_sample_treat)
  } else {
    out$train_predictions_control <- avg("train_predictions_control")
    out$test_predictions_control  <- avg("test_predictions_control")
    out$train_predictions_treat   <- avg("train_predictions_treat")
    out$test_predictions_treat    <- avg("test_predictions_treat")
  }
  # Total is always averaged per-chain (control + treat are on different scales)
  out$train_predictions <- avg("train_predictions")
  out$test_predictions  <- avg("test_predictions")

  # ── Sigma ─────────────────────────────────────────────────────────────────
  out$sigma <- pool_vector("sigma")

  # ── Acceptance ratios ─────────────────────────────────────────────────────
  out$acceptance_ratio_control <- pool_vector("acceptance_ratio_control")
  out$acceptance_ratio_treat   <- pool_vector("acceptance_ratio_treat")

  # ── Variable importance ───────────────────────────────────────────────────
  out$split_probs_control        <- pool_matrix("split_probs_control")
  out$store_split_counts_control <- pool_matrix("store_split_counts_control")
  out$split_probs_treat          <- pool_matrix("split_probs_treat")
  out$store_split_counts_treat   <- pool_matrix("store_split_counts_treat")

  # ── Chain-level diagnostics ───────────────────────────────────────────────
  out$chains <- list(
    n_chains                  = n_chains,
    acceptance_ratios_control = vapply(chains,
      function(ch) mean(ch$acceptance_ratio_control), numeric(1)),
    acceptance_ratios_treat   = vapply(chains,
      function(ch) mean(ch$acceptance_ratio_treat),   numeric(1))
  )

  out$mcmc$n_chains <- n_chains
  out
}

#' Posterior ATE samples from a CATE draw matrix
#'
#' @param tau_samples S x n matrix of posterior CATE draws.
#' @param bayesian_bootstrap If \code{TRUE}, reweight columns at each iteration
#'   with Dirichlet(1, ..., 1) weights, giving posterior draws of the
#'   population ATE. If \code{FALSE}, use equal 1/n weights (mixed ATE,
#'   conditional on the observed covariates).
#' @return Length-S numeric vector of posterior ATE draws.
#' @importFrom stats rexp
#' @noRd
.ate_samples <- function(tau_samples, bayesian_bootstrap = TRUE) {
  if (!bayesian_bootstrap) return(rowMeans(tau_samples))
  S <- nrow(tau_samples)
  n <- ncol(tau_samples)
  W <- matrix(rexp(S * n), nrow = S, ncol = n)
  W <- W / rowSums(W)
  rowSums(W * tau_samples)
}


# Error message for the removed forest-wide horseshoe prior ------------------
#
# The forest-wide update was never run by the sampler, so "horseshoe_fw" fitted
# the same model as "half-cauchy" with the same local scale. The option was
# removed in 2.1.1; this message points users to the equivalent call.
.horseshoe_fw_removed <- function(arg) {
  local_arg <- sub("prior_type", "local_hp", arg, fixed = TRUE)
  paste0(arg, " = 'horseshoe_fw' is no longer supported. ",
         "Use ", arg, " = 'half-cauchy' with the same ", local_arg,
         "; it fits the same model.")
}


# Calibrated default horseshoe scale, single-forest models -------------------
#
# Used by HorseTrees() and ShrinkageTrees(). The local and global scales are
# both set to k / sqrt(number_of_trees).
#
# The leaf standard deviation is sqrt(omega) * local_hp * global_hp, and a sum
# of m heavy-tailed leaves grows like m, so the implied prior scale of the
# whole ensemble is sqrt(omega) * k^2 -- free of number_of_trees. The value
# below puts that quantity near 1 on a standardized response, which gave
# near-nominal pointwise coverage for p between 50 and 5000.
#
# Keep in sync with the `k` default in the HorseTrees() signature; the
# regression test in tests/testthat/test-horseshoe-defaults.R asserts this.
HORSESHOE_K_SINGLE <- 1.0

# Calibrated default horseshoe scale, causal (tau-learner) models ------------
#
# Used by CausalHorseForest() and CausalShrinkageForest(), applied to both
# forests using each forest's own tree count.
#
# Larger than HORSESHOE_K_SINGLE because the treatment forest enters as
# b * tau(x) with b = +/- 1/2, so its contribution to the response is halved.
# Measured in terms of the contribution scale the two defaults agree closely.
#
# Keep in sync with the `k` default in the CausalHorseForest() signature.
HORSESHOE_K_CAUSAL <- 1.5
