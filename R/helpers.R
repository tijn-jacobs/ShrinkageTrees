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
#' @param ic_indicator Optional numeric vector; 1 = interval-censored, 0 = right-censored/exact.
#'
#' @return A list with elements:
#'   \item{mu}{Estimated mean of survival times.}
#'   \item{sd}{Estimated standard deviation of survival times.}
#'   \item{min}{Estimated minimum of survival times.}
#'   \item{max}{Estimated maximum of survival times.}
#'
#' @importFrom stats dnorm pnorm
censored_info <- function(y, status, left_time = NULL, right_time = NULL,
                          ic_indicator = NULL) {

  if (requireNamespace("survival", quietly = TRUE)) {

    if (!is.null(left_time) && !is.null(right_time)) {
      # Interval-censored path: use Surv(time, time2, type = "interval2")
      surv_obj <- survival::Surv(left_time, right_time, type = "interval2")
      fit <- survival::survreg(surv_obj ~ 1, dist = "gaussian")

      mu <- as.numeric(fit$coefficients[1])
      sd <- as.numeric(fit$scale)

      # Impute censored values for min/max estimation
      imputed_y <- y

      # Right-censored observations: impute with Mills ratio
      if (!is.null(ic_indicator)) {
        rc_idx <- which(status == 0 & ic_indicator == 0)
        if (length(rc_idx) > 0) {
          a <- (y[rc_idx] - mu) / sd
          lambda_val <- dnorm(a) / (1 - pnorm(a))
          imputed_y[rc_idx] <- mu + sd * lambda_val
        }
      }

      min_val <- min(imputed_y, na.rm = TRUE)
      max_val <- max(imputed_y, na.rm = TRUE)

    } else {
      # Right-censored path (existing)
      fit <- survival::survreg(survival::Surv(y, status) ~ 1, dist = "gaussian")

      mu <- as.numeric(fit$coefficients[1])
      sd <- as.numeric(fit$scale)

      # ---- IMPUTE censored values ----
      cens_idx <- which(status == 0)
      imputed_y <- y

      if (length(cens_idx) > 0) {
        a <- (y[cens_idx] - mu) / sd
        lambda <- dnorm(a) / (1 - pnorm(a))     # Mills ratio
        imputed_y[cens_idx] <- mu + sd * lambda # Replace censored with E[Y | Y > c]
      }

      # min/max after imputation
      min_val <- min(imputed_y, na.rm = TRUE)
      max_val <- max(imputed_y, na.rm = TRUE)
    }

  } else {

    # Fallback: use only uncensored observations
    min_val <- min(y, na.rm = TRUE)
    max_val <- max(y, na.rm = TRUE)

    mu <- mean(y[status == 1], na.rm = TRUE)
    sd <- sd(y[status == 1], na.rm = TRUE)
  }

  return(list(
    mu = mu,
    sd = sd,
    min = min_val,
    max = max_val
  ))
}

# Internal helper: pool a list of ShrinkageTrees objects (one per chain) into
# a single ShrinkageTrees object with concatenated/averaged posteriors.
.combine_chains <- function(chains) {
  n_chains <- length(chains)
  out      <- chains[[1]]   # use first chain as structural template

  pool_matrix <- function(field) {
    mats <- lapply(chains, `[[`, field)
    if (all(vapply(mats, is.null, logical(1)))) return(NULL)
    do.call(rbind, mats)
  }
  pool_vector <- function(field) {
    vecs <- lapply(chains, `[[`, field)
    if (all(vapply(vecs, is.null, logical(1)))) return(NULL)
    unlist(vecs)
  }

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

  pool_matrix <- function(field) {
    mats <- lapply(chains, `[[`, field)
    if (all(vapply(mats, is.null, logical(1)))) return(NULL)
    do.call(rbind, mats)
  }
  pool_vector <- function(field) {
    vecs <- lapply(chains, `[[`, field)
    if (all(vapply(vecs, is.null, logical(1)))) return(NULL)
    unlist(vecs)
  }
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
