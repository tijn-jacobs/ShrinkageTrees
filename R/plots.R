# ── Internal helpers ──────────────────────────────────────────────────────────

.check_bayesplot <- function() {
  if (!requireNamespace("bayesplot", quietly = TRUE)) {
    stop(
      "Package 'bayesplot' is needed for plotting. ",
      "Install it with: install.packages('bayesplot')",
      call. = FALSE
    )
  }
}

.check_ggplot2 <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "Package 'ggplot2' is needed for this plot type. ",
      "Install it with: install.packages('ggplot2')",
      call. = FALSE
    )
  }
}

# Reshape a concatenated posterior vector into an [iterations x chains x 1]
# array accepted by bayesplot trace/density functions.
.to_draws_array <- function(vec, par_name, N_post, n_chains) {
  array(
    matrix(vec, nrow = N_post, ncol = n_chains),
    dim      = c(N_post, n_chains, 1L),
    dimnames = list(
      NULL,
      paste0("chain:", seq_len(n_chains)),
      par_name
    )
  )
}

# Extract feature names from a matrix, falling back to X1, X2, ...
.feature_names <- function(mat) {
  nms <- colnames(mat)
  if (is.null(nms)) nms <- paste0("X", seq_len(ncol(mat)))
  nms
}

# Build a top-n variable importance plot from a split_probs matrix.
.vi_intervals <- function(sp, X_mat, n_vi, title, ...) {
  colnames(sp) <- .feature_names(X_mat)
  top <- order(colMeans(sp), decreasing = TRUE)[seq_len(min(n_vi, ncol(sp)))]
  bayesplot::mcmc_intervals(
    sp[, top, drop = FALSE],
    prob       = 0.5,
    prob_outer = 0.95,
    ...
  ) + ggplot2::ggtitle(title)
}


# ── plot.ShrinkageTrees ───────────────────────────────────────────────────────

#' Plot diagnostics for a ShrinkageTrees model
#'
#' Visualises posterior draws using \pkg{bayesplot}. Requires the suggested
#' packages \pkg{bayesplot} and \pkg{ggplot2}.
#'
#' @param x A \code{ShrinkageTrees} object.
#' @param type Character; one of:
#'   \describe{
#'     \item{\code{"trace"}}{Sigma traceplot across MCMC iterations (one line
#'       per chain). Useful for assessing chain mixing.}
#'     \item{\code{"density"}}{Overlaid posterior density of sigma, one curve
#'       per chain.}
#'     \item{\code{"vi"}}{Posterior credible intervals for variable inclusion
#'       probabilities (top \code{n_vi} predictors). Only available for
#'       Dirichlet prior models.}
#'   }
#' @param n_vi Integer; number of top variables to display when
#'   \code{type = "vi"}. Default \code{10}.
#' @param ... Additional arguments forwarded to the underlying
#'   \pkg{bayesplot} function.
#' @return A \pkg{ggplot2} object.
#'
#' @examples
#' \dontrun{
#' # Fit a model with multiple chains
#' fit <- ShrinkageTrees(
#'   y = y, X_train = X,
#'   prior_type = "horseshoe",
#'   local_hp = 0.1, global_hp = 0.1,
#'   N_post = 500, N_burn = 200,
#'   n_chains = 4
#' )
#'
#' # Sigma traceplot — check chain mixing
#' plot(fit, type = "trace")
#'
#' # Overlaid posterior densities of sigma per chain
#' plot(fit, type = "density")
#'
#' # Variable importance (Dirichlet prior only)
#' plot(fit, type = "vi", n_vi = 10)
#' }
#' @export
plot.ShrinkageTrees <- function(x,
                                type = c("trace", "density", "vi"),
                                n_vi = 10,
                                ...) {
  type <- match.arg(type)
  .check_bayesplot()
  .check_ggplot2()

  n_ch   <- if (!is.null(x$mcmc$n_chains)) x$mcmc$n_chains else 1L
  N_post <- x$mcmc$N_post

  if (type == "trace") {
    arr <- .to_draws_array(x$sigma, "sigma", N_post, n_ch)
    bayesplot::mcmc_trace(arr, ...)

  } else if (type == "density") {
    arr <- .to_draws_array(x$sigma, "sigma", N_post, n_ch)
    bayesplot::mcmc_dens_overlay(arr, ...)

  } else {
    sp <- x$split_probs
    if (is.null(sp))
      stop("Variable importance not available (split_probs is NULL).",
           call. = FALSE)
    .vi_intervals(sp, x$data$X_train, n_vi,
                  title = "Variable importance", ...)
  }
}


# ── plot.CausalShrinkageForest ────────────────────────────────────────────────

#' Plot diagnostics for a CausalShrinkageForest model
#'
#' Visualises posterior draws using \pkg{bayesplot} and \pkg{ggplot2}.
#' Requires the suggested packages \pkg{bayesplot} and \pkg{ggplot2}.
#'
#' @param x A \code{CausalShrinkageForest} object.
#' @param type Character; one of:
#'   \describe{
#'     \item{\code{"trace"}}{Sigma traceplot (chain mixing).}
#'     \item{\code{"density"}}{Overlaid posterior density of sigma per chain.}
#'     \item{\code{"ate"}}{Posterior density of the average treatment effect
#'       (ATE) with 95 \% credible region. Requires
#'       \code{store_posterior_sample = TRUE}.}
#'     \item{\code{"cate"}}{Point estimates and 95 \% credible intervals for
#'       the CATE of each training observation, sorted by posterior mean.
#'       Requires \code{store_posterior_sample = TRUE}.}
#'     \item{\code{"vi"}}{Posterior credible intervals for variable inclusion
#'       probabilities (Dirichlet prior only). Controlled by \code{forest}.}
#'   }
#' @param forest For \code{type = "vi"}: which forest to display.
#'   One of \code{"both"} (default), \code{"control"}, or \code{"treat"}.
#'   When \code{"both"}, a named list of two \pkg{ggplot2} objects is returned.
#' @param n_vi Integer; number of top variables for \code{type = "vi"}.
#'   Default \code{10}.
#' @param ... Additional arguments forwarded to the underlying
#'   \pkg{bayesplot} or \pkg{ggplot2} function.
#' @return A \pkg{ggplot2} object, or (for \code{type = "vi"} with
#'   \code{forest = "both"}) a named list with elements \code{control} and
#'   \code{treat}.
#'
#' @examples
#' \dontrun{
#' # Fit a causal forest with multiple chains
#' fit <- CausalShrinkageForest(
#'   y = y,
#'   X_train_control = X, X_train_treat = X,
#'   treatment_indicator_train = w,
#'   prior_type_control = "horseshoe", prior_type_treat = "horseshoe",
#'   local_hp_control = 0.1, global_hp_control = 0.1,
#'   local_hp_treat  = 0.1, global_hp_treat  = 0.1,
#'   N_post = 500, N_burn = 200,
#'   store_posterior_sample = TRUE,
#'   n_chains = 4
#' )
#'
#' # Sigma traceplot — check chain mixing
#' plot(fit, type = "trace")
#'
#' # Overlaid posterior densities of sigma
#' plot(fit, type = "density")
#'
#' # Posterior distribution of the average treatment effect
#' plot(fit, type = "ate")
#'
#' # Conditional treatment effects for each observation
#' plot(fit, type = "cate")
#'
#' # Variable importance for both forests (Dirichlet prior only)
#' vi_plots <- plot(fit, type = "vi", forest = "both")
#' vi_plots$control
#' vi_plots$treat
#' }
#' @export
plot.CausalShrinkageForest <- function(x,
                                       type   = c("trace", "density",
                                                  "ate", "cate", "vi"),
                                       forest = c("both", "control", "treat"),
                                       n_vi   = 10,
                                       ...) {
  type   <- match.arg(type)
  forest <- match.arg(forest)
  .check_bayesplot()
  .check_ggplot2()

  n_ch   <- if (!is.null(x$mcmc$n_chains)) x$mcmc$n_chains else 1L
  N_post <- x$mcmc$N_post

  # ── Sigma trace / density ────────────────────────────────────────────────
  if (type == "trace") {
    arr <- .to_draws_array(x$sigma, "sigma", N_post, n_ch)
    return(bayesplot::mcmc_trace(arr, ...))
  }

  if (type == "density") {
    arr <- .to_draws_array(x$sigma, "sigma", N_post, n_ch)
    return(bayesplot::mcmc_dens_overlay(arr, ...))
  }

  # ── ATE posterior ────────────────────────────────────────────────────────
  if (type == "ate") {
    sc <- x$train_predictions_sample_control
    st <- x$train_predictions_sample_treat
    if (is.null(sc) || is.null(st))
      stop("ATE plot requires store_posterior_sample = TRUE.", call. = FALSE)

    ate_draws <- rowMeans(st) - rowMeans(sc)
    arr <- .to_draws_array(ate_draws, "ATE", N_post, n_ch)
    return(bayesplot::mcmc_areas(arr, prob = 0.95, ...))
  }

  # ── CATE distribution ────────────────────────────────────────────────────
  if (type == "cate") {
    sc <- x$train_predictions_sample_control
    st <- x$train_predictions_sample_treat
    if (is.null(sc) || is.null(st))
      stop("CATE plot requires store_posterior_sample = TRUE.", call. = FALSE)

    cate_draws <- st - sc
    cate_mean  <- colMeans(cate_draws)
    cate_lo    <- apply(cate_draws, 2, quantile, 0.025)
    cate_hi    <- apply(cate_draws, 2, quantile, 0.975)

    ord <- order(cate_mean)
    df  <- data.frame(
      obs  = seq_along(cate_mean),
      cate = cate_mean[ord],
      lo   = cate_lo[ord],
      hi   = cate_hi[ord]
    )

    return(
      ggplot2::ggplot(df, ggplot2::aes(x = .data$cate, y = .data$obs)) +
        ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                            colour = "grey50") +
        ggplot2::geom_linerange(
          ggplot2::aes(xmin = .data$lo, xmax = .data$hi),
          colour = "grey70"
        ) +
        ggplot2::geom_point(size = 0.8) +
        ggplot2::labs(
          x     = "CATE",
          y     = "Observation (sorted by posterior mean CATE)",
          title = "Conditional Average Treatment Effects"
        ) +
        bayesplot::theme_default()
    )
  }

  # ── Variable importance ──────────────────────────────────────────────────
  if (forest %in% c("control", "both")) {
    sp_c <- x$split_probs_control
    if (is.null(sp_c))
      stop("Control variable importance not available.", call. = FALSE)
    p_control <- .vi_intervals(sp_c, x$data$X_train_control, n_vi,
                               title = "Variable importance — control forest",
                               ...)
  }
  if (forest %in% c("treat", "both")) {
    sp_t <- x$split_probs_treat
    if (is.null(sp_t))
      stop("Treatment variable importance not available.", call. = FALSE)
    p_treat <- .vi_intervals(sp_t, x$data$X_train_treat, n_vi,
                             title = "Variable importance — treatment forest",
                             ...)
  }

  if (forest == "control") return(p_control)
  if (forest == "treat")   return(p_treat)
  list(control = p_control, treat = p_treat)
}
