#' Horseshoe Regression Trees (HorseTrees)
#'
#' Fits a Bayesian Horseshoe Trees model with a single learner.
#' Implements regularization on the step heights using a global-local Horseshoe
#' prior, controlled via the parameter \code{k}. Supports continuous, binary,
#' right-censored, and interval-censored (survival) outcomes.
#'
#' For continuous outcomes, the model centers and optionally standardizes the
#' outcome using a prior guess of the standard deviation.
#' For binary outcomes, the function uses a probit link formulation.
#' For right-censored outcomes (survival data), the function can handle
#' follow-up times either on the original time scale or log-transformed.
#' For interval-censored outcomes, provide \code{left_time} and
#' \code{right_time} instead of \code{y} and \code{status}; the event
#' indicators are derived internally following the
#' \code{survival::Surv(type = "interval2")} convention.
#' Generalized implementation with multiple prior possibilities is given by
#' \code{\link{ShrinkageTrees}}.
#' 
#' @examples
#' # Minimal example: continuous outcome
#' n <- 25
#' p <- 5
#' X <- matrix(rnorm(n * p), ncol = p)
#' y <- X[, 1] + rnorm(n)
#' fit1 <- HorseTrees(y = y, X_train = X, outcome_type = "continuous", 
#'                    number_of_trees = 5, N_post = 75, N_burn = 25, 
#'                    verbose = FALSE)
#'
#' # Minimal example: binary outcome
#' X <- matrix(rnorm(n * p), ncol = p)
#' y <- ifelse(X[, 1] + rnorm(n) > 0, 1, 0)
#' fit2 <- HorseTrees(y = y, X_train = X, outcome_type = "binary", 
#'                    number_of_trees = 5, N_post = 75, N_burn = 25, 
#'                    verbose = FALSE)
#'
#' # Minimal example: right-censored outcome
#' X <- matrix(rnorm(n * p), ncol = p)
#' time <- rexp(n, rate = 0.1)
#' status <- rbinom(n, 1, 0.7)
#' fit3 <- HorseTrees(y = time, status = status, X_train = X,
#'                    outcome_type = "right-censored", number_of_trees = 5,
#'                    N_post = 75, N_burn = 25, verbose = FALSE)
#'
#' # Minimal example: interval-censored outcome
#' X <- matrix(rnorm(n * p), ncol = p)
#' true_t <- rexp(n, rate = 0.1)
#' left_t  <- true_t * runif(n, 0.5, 1)
#' right_t <- true_t * runif(n, 1, 1.5)
#' # Mark some as exact, some as right-censored
#' exact <- sample(n, 8); left_t[exact] <- true_t[exact]; right_t[exact] <- true_t[exact]
#' rc <- sample(setdiff(seq_len(n), exact), 5); right_t[rc] <- Inf
#' fit4 <- HorseTrees(left_time = left_t, right_time = right_t, X_train = X,
#'                    outcome_type = "interval-censored", number_of_trees = 5,
#'                    N_post = 75, N_burn = 25, verbose = FALSE)
#'
#' # Larger continuous example (not run automatically)
#' \donttest{
#' n <- 100
#' p <- 100
#' X <- matrix(rnorm(100 * p), ncol = p)
#' X_test <- matrix(rnorm(50 * p), ncol = p)
#' y <- X[, 1] + X[, 2] - X[, 3] + rnorm(100, sd = 0.5)
#'
#' fit5 <- HorseTrees(y = y,
#'                    X_train = X,
#'                    X_test = X_test,
#'                    outcome_type = "continuous",
#'                    number_of_trees = 200,
#'                    N_post = 2500,
#'                    N_burn = 2500,
#'                    store_posterior_sample = TRUE,
#'                    verbose = TRUE)
#'
#' plot(fit4$sigma, type = "l", ylab = expression(sigma),
#'      xlab = "Iteration", main = "Sigma traceplot")
#'
#' hist(fit4$train_predictions_sample[, 1],
#'      main = "Posterior distribution of prediction outcome individual 1",
#'      xlab = "Prediction", breaks = 20)
#' }
#'                        
#' @param y Outcome vector. Numeric. Can represent continuous outcomes, binary
#' outcomes (0/1), or follow-up times for survival data. Set to \code{NULL}
#' (default) when using \code{outcome_type = "interval-censored"}, as values
#' are derived from \code{left_time} and \code{right_time}.
#' @param status Optional censoring indicator vector (1 = event occurred,
#' 0 = censored). Required if \code{outcome_type = "right-censored"}.
#' For interval-censored outcomes, this is derived automatically from
#' \code{left_time} and \code{right_time}.
#' @param X_train Covariate matrix for training. Each row corresponds to an
#' observation, and each column to a covariate.
#' @param X_test Optional covariate matrix for test data. If NULL, defaults to
#' the mean of the training covariates.
#' @param left_time Optional numeric vector of left (lower) time boundaries.
#' Required when \code{outcome_type = "interval-censored"}. Exact events
#' have \code{left_time == right_time}; right-censored observations have
#' \code{right_time = Inf}; interval-censored observations have finite
#' \code{left_time < right_time}.
#' @param right_time Optional numeric vector of right (upper) time boundaries.
#' Required when \code{outcome_type = "interval-censored"}. Use \code{Inf}
#' for right-censored observations.
#' @param outcome_type Type of outcome. One of \code{"continuous"},
#' \code{"binary"}, \code{"right-censored"}, or \code{"interval-censored"}.
#' @param timescale Indicates the scale of follow-up times. Options are 
#' `"time"` (nonnegative follow-up times, will be log-transformed internally) 
#'   or `"log"` (already log-transformed). Only used when 
#'   `outcome_type = "right-censored"` or `"interval-censored"`.
#' @param number_of_trees Number of trees in the ensemble. Default is 200.
#' @param k Horseshoe scale hyperparameter. Default is 1.0. Controls the
#' overall level of shrinkage by setting the scale of both the local and the
#' global component, parameterized as
#' \eqn{\alpha_{\mathrm{local}} = \alpha_{\mathrm{global}} =
#' k / \sqrt{\mathrm{number\_of\_trees}}}.
#'
#' \strong{Choosing k.} Values between roughly 0.5 and 1.5 worked well across a
#' wide range of simulated settings. Within that range:
#' \itemize{
#'   \item \strong{smaller k shrinks more aggressively} -- tighter credible
#'     intervals, more leaves pulled to zero. Prefer it when the signal is
#'     sparse, the dimension is high, or the sample is small.
#'   \item \strong{larger k is more conservative} -- wider intervals, more
#'     freedom for individual leaves. Prefer it when the function is complex
#'     or the signal is strong.
#' }
#' The default sits in the middle. If you would rather select it from the data,
#' \code{k} is a natural target for cross-validation; see the "Choosing k"
#' section of the package vignette.
#'
#' \strong{Changed in version 2.1.0.} The default was 0.1 in earlier releases.
#' That value was calibrated against a defect in the global update which made
#' \code{global_hp} inert (see \code{NEWS.md}); with the defect corrected it
#' shrinks far too aggressively. Pass \code{k = 0.1} explicitly only if you
#' need to reproduce pre-2.1.0 output, and note that this does not reproduce it
#' exactly, since the sampler itself has been corrected.
#' @param power Power parameter for tree structure prior. Default is 2.0.
#' @param base Base parameter for tree structure prior. Default is 0.95.
#' @param p_grow Probability of proposing a grow move. Default is 0.4.
#' @param p_prune Probability of proposing a prune move. Default is 0.4.
#' @param nu Degrees of freedom for the error distribution prior. Default is 3.
#' @param q Quantile hyperparameter for the error variance prior. 
#' Default is 0.90.
#' @param sigma Optional known value for error standard deviation. If NULL, 
#' estimated from data.
#' @param N_post Number of posterior samples to store. Default is 1000.
#' @param N_burn Number of burn-in iterations. Default is 1000.
#' @param delayed_proposal Number of delayed iterations before proposal. Only 
#' for reversible updates. Default is 5.
#' @param store_posterior_sample Logical; whether to store posterior samples for
#'  each iteration. Default is TRUE.
#' @param n_chains Number of independent MCMC chains to run. Default is
#'   \code{1} (standard single-chain behaviour). When \code{n_chains > 1} the
#'   chains are run in parallel via \code{parallel::mclapply} and their
#'   posterior samples are pooled into a single \code{ShrinkageTrees} object,
#'   so all existing \code{print}, \code{summary}, and \code{predict} methods
#'   work without modification. On Windows, \code{mclapply} falls back to
#'   sequential execution.
#' @param verbose Logical; whether to print verbose output. Default is TRUE.
#'
#' @return An S3 object of class \code{"ShrinkageTrees"} with the following elements:
#' \describe{
#'   \item{train_predictions}{Vector of posterior mean predictions on the 
#'   training data.}
#'   \item{test_predictions}{Vector of posterior mean predictions on the test 
#'   data (or on mean covariate vector if `X_test` not provided).}
#'   \item{sigma}{Vector of posterior samples of the error variance.}
#'   \item{acceptance_ratio}{Average acceptance ratio across trees during 
#'   sampling.}
#'   \item{train_predictions_sample}{Matrix of posterior samples of training 
#'   predictions (iterations in rows, observations in columns). Present only 
#'   if `store_posterior_sample = TRUE`.}
#'   \item{test_predictions_sample}{Matrix of posterior samples of test 
#'   predictions. Present only if `store_posterior_sample = TRUE`.}
#'   \item{train_probabilities}{Vector of posterior mean probabilities on the 
#'   training data (only for `outcome_type = "binary"`).}
#'   \item{test_probabilities}{Vector of posterior mean probabilities on the 
#'   test data (only for `outcome_type = "binary"`).}
#'   \item{train_probabilities_sample}{Matrix of posterior samples of training 
#'   probabilities (only for `outcome_type = "binary"` and if 
#'   `store_posterior_sample = TRUE`).}
#'   \item{test_probabilities_sample}{Matrix of posterior samples of test 
#'   probabilities (only for `outcome_type = "binary"` and if 
#'   `store_posterior_sample = TRUE`).}
#' }
#' 
#' @seealso
#' Model family: \code{\link{ShrinkageTrees}} (flexible prior choice),
#' \code{\link{CausalHorseForest}} (causal inference),
#' \code{\link{CausalShrinkageForest}} (causal, flexible prior).
#'
#' Survival wrappers: \code{\link{SurvivalBART}}, \code{\link{SurvivalDART}}.
#'
#' S3 methods: \code{\link{print.ShrinkageTrees}},
#' \code{\link{summary.ShrinkageTrees}},
#' \code{\link{predict.ShrinkageTrees}},
#' \code{\link{plot.ShrinkageTrees}}.
#' 
#' @importFrom Rcpp evalCpp
#' @useDynLib ShrinkageTrees, .registration = TRUE
#' @importFrom stats sd qchisq qnorm pnorm
#' @importFrom parallel mclapply detectCores
#' @export

HorseTrees <- function(y = NULL,
                       status = NULL,
                       X_train,
                       X_test = NULL,
                       left_time = NULL,
                       right_time = NULL,
                       outcome_type = "continuous",
                       timescale = "time",
                       number_of_trees = 200,
                       k = 1.0,
                       power = 2.0,
                       base = 0.95,
                       p_grow = 0.4,
                       p_prune = 0.4,
                       nu = 3,
                       q = 0.90,
                       sigma = NULL,
                       N_post = 1000,
                       N_burn = 1000,
                       delayed_proposal = 5,
                       store_posterior_sample = TRUE,
                       n_chains = 1,
                       verbose = TRUE) {

  # HorseTrees() is ShrinkageTrees() with the horseshoe prior and both scales
  # set to k / sqrt(number_of_trees).
  hp <- k / sqrt(number_of_trees)

  fit <- ShrinkageTrees(
    y                      = y,
    status                 = status,
    X_train                = X_train,
    X_test                 = X_test,
    left_time              = left_time,
    right_time             = right_time,
    outcome_type           = outcome_type,
    timescale              = timescale,
    number_of_trees        = number_of_trees,
    prior_type             = "horseshoe",
    local_hp               = hp,
    global_hp              = hp,
    power                  = power,
    base                   = base,
    p_grow                 = p_grow,
    p_prune                = p_prune,
    nu                     = nu,
    q                      = q,
    sigma                  = sigma,
    N_post                 = N_post,
    N_burn                 = N_burn,
    delayed_proposal       = delayed_proposal,
    store_posterior_sample = store_posterior_sample,
    n_chains               = n_chains,
    verbose                = verbose
  )

  fit$call <- match.call()
  fit
}
