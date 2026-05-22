#' Bayesian bootstrap average treatment effect
#'
#' Post-hoc reweights the stored posterior CATE draws of a fitted causal
#' model to produce credible intervals for the \emph{population} ATE (PATE)
#' that incorporate uncertainty in the covariate distribution
#' \eqn{F_X}.
#'
#' At each MCMC iteration \eqn{s} the conditional treatment effects
#' \eqn{\tau^{(s)}(x_i)} are reweighted with
#' \eqn{(w_1^{(s)}, \dots, w_n^{(s)}) \sim \mathrm{Dir}(1, \dots, 1)}
#' to give a draw
#' \deqn{\widehat{\mathrm{PATE}}^{(s)} = \sum_{i=1}^n w_i^{(s)}\, \tau^{(s)}(x_i).}
#' The collection \eqn{\{\widehat{\mathrm{PATE}}^{(s)}\}} approximates the
#' posterior of the PATE, integrating over \eqn{\tau(\cdot)} \emph{and}
#' \eqn{F_X}. The equal-weight mixed ATE (MATE),
#' \eqn{\widehat{\mathrm{MATE}}^{(s)} = n^{-1}\sum_i \tau^{(s)}(x_i)},
#' is returned alongside for comparison.
#'
#' For reproducibility, call \code{set.seed()} before invoking the function
#' to fix the Dirichlet draws.
#'
#' @param object Either a fitted \code{CausalShrinkageForest} (from
#'   \code{\link{CausalShrinkageForest}} or \code{\link{CausalHorseForest}}
#'   with \code{store_posterior_sample = TRUE}) or a
#'   \code{CausalShrinkageForestPrediction} (from
#'   \code{\link{predict.CausalShrinkageForest}}).
#' @param alpha One minus the credible level. Default \code{0.05}
#'   (a 95 percent credible interval).
#' @return A list with
#'   \describe{
#'     \item{pate_mean, pate_ci, pate_samples}{Posterior mean, credible
#'       interval (named \code{lower} and \code{upper}), and full vector of
#'       draws of the Bayesian-bootstrap PATE.}
#'     \item{mate_mean, mate_ci, mate_samples}{Same quantities for the
#'       equal-weight mixed ATE.}
#'     \item{n, S}{Number of observations and posterior draws used.}
#'   }
#' @seealso \code{\link{summary.CausalShrinkageForest}},
#'   \code{\link{plot.CausalShrinkageForest}},
#'   \code{\link{predict.CausalShrinkageForest}}
#'
#' @examples
#' # Small toy causal model (binary outcome, for speed)
#' set.seed(1)
#' n <- 40; p <- 3
#' X <- matrix(runif(n * p), ncol = p)
#' trt <- rbinom(n, 1, 0.5)
#' y <- X[, 1] + trt * (0.5 + X[, 2]) + rnorm(n)
#'
#' fit <- CausalShrinkageForest(
#'   y = y,
#'   X_train_control = X, X_train_treat = X,
#'   treatment_indicator_train = trt,
#'   outcome_type = "continuous",
#'   number_of_trees_control = 5, number_of_trees_treat = 5,
#'   prior_type_control = "horseshoe", prior_type_treat = "horseshoe",
#'   local_hp_control = 0.1, global_hp_control = 0.1,
#'   local_hp_treat = 0.1, global_hp_treat = 0.1,
#'   N_post = 20, N_burn = 10,
#'   store_posterior_sample = TRUE,
#'   verbose = FALSE
#' )
#'
#' bb <- bayesian_bootstrap_ate(fit, alpha = 0.05)
#' bb$pate_mean
#' bb$pate_ci
#'
#' @export
bayesian_bootstrap_ate <- function(object, alpha = 0.05) {

  if (inherits(object, "CausalShrinkageForest")) {
    tau_samples <- object$train_predictions_sample_treat
    if (is.null(tau_samples))
      stop("Posterior CATE samples not stored. Refit with ",
           "store_posterior_sample = TRUE.", call. = FALSE)
  } else if (inherits(object, "CausalShrinkageForestPrediction")) {
    tau_samples <- object$cate_samples
    if (is.null(tau_samples))
      stop("Prediction object does not contain cate_samples. ",
           "Refit and re-predict with the current package version.",
           call. = FALSE)
  } else {
    stop("`object` must be a CausalShrinkageForest or ",
         "CausalShrinkageForestPrediction.", call. = FALSE)
  }

  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0 || alpha >= 1)
    stop("`alpha` must be a single number in (0, 1).", call. = FALSE)

  probs <- c(alpha / 2, 1 - alpha / 2)

  mate_samples <- .ate_samples(tau_samples, bayesian_bootstrap = FALSE)
  pate_samples <- .ate_samples(tau_samples, bayesian_bootstrap = TRUE)

  ci <- function(s) {
    q <- quantile(s, probs = probs)
    list(lower = unname(q[1]), upper = unname(q[2]))
  }

  list(
    pate_mean    = mean(pate_samples),
    pate_ci      = ci(pate_samples),
    pate_samples = pate_samples,
    mate_mean    = mean(mate_samples),
    mate_ci      = ci(mate_samples),
    mate_samples = mate_samples,
    n            = ncol(tau_samples),
    S            = nrow(tau_samples)
  )
}
