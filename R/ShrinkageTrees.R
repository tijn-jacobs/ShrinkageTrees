#' General Shrinkage Regression Trees (ShrinkageTrees)
#'
#' Fits a Bayesian Shrinkage Tree model with flexible global-local priors on the
#' step heights. This function generalizes \code{\link{HorseTrees}} by allowing
#' different global-local shrinkage priors on the step heights.
#'
#' @param y Outcome vector. Numeric. Can represent continuous outcomes, binary
#' outcomes (0/1), or follow-up times for survival data.
#' @param status Optional censoring indicator vector (1 = event occurred,
#' 0 = censored). Required if `outcome_type = "right-censored"`.
#' @param X_train Covariate matrix for training. Each row corresponds to an
#' observation, and each column to a covariate.
#' @param X_test Optional covariate matrix for test data. If NULL, defaults to
#' the mean of the training covariates.
#' @param outcome_type Type of outcome. One of `"continuous"`, `"binary"`, or
#' `"right-censored"`.
#' @param timescale Indicates the scale of follow-up times. Options are
#' `"time"` (nonnegative follow-up times, will be log-transformed internally)
#' or `"log"` (already log-transformed). Only used when
#' `outcome_type = "right-censored"`.
#' @param number_of_trees Number of trees in the ensemble. Default is 200.
#' @param prior_type Type of prior on the step heights. Options include
#' `"horseshoe"`, `"horseshoe_fw"`, `"horseshoe_EB"`, and `"half-cauchy"`.
#' @param local_hp Local hyperparameter controlling shrinkage on individual
#' step heights. Should typically be set smaller than 1 / sqrt(number_of_trees).
#' @param global_hp Global hyperparameter controlling overall shrinkage.
#' Must be specified for Horseshoe-type priors; ignored for `prior_type = "half-cauchy"`.
#' @param power Power parameter for the tree structure prior. Default is 2.0.
#' @param base Base parameter for the tree structure prior. Default is 0.95.
#' @param p_grow Probability of proposing a grow move. Default is 0.4.
#' @param p_prune Probability of proposing a prune move. Default is 0.4.
#' @param nu Degrees of freedom for the error distribution prior. Default is 3.
#' @param q Quantile hyperparameter for the error variance prior. Default is 0.90.
#' @param sigma Optional known value for error standard deviation. If NULL,
#' estimated from data.
#' @param N_post Number of posterior samples to store. Default is 1000.
#' @param N_burn Number of burn-in iterations. Default is 1000.
#' @param delayed_proposal Number of delayed iterations before proposal. Only
#' for reversible updates. Default is 5.
#' @param store_posterior_sample Logical; whether to store posterior samples for
#' each iteration. Default is TRUE.
#' @param seed Random seed for reproducibility.
#' @param verbose Logical; whether to print verbose output. Default is TRUE.
#'
#' @return A named list with the following elements:
#' \describe{
#'   \item{train_predictions}{Vector of posterior mean predictions on the
#'   training data.}
#'   \item{test_predictions}{Vector of posterior mean predictions on the test
#'   data (or on mean covariate vector if `X_test` not provided).}
#'   \item{sigma}{Vector of posterior samples of the error variance.}
#'   \item{acceptance_ratio}{Average acceptance ratio across trees during
#'   sampling.}
#'   \item{train_predictions_sample}{Matrix of posterior samples of training
#'   predictions (iterations in rows, observations in columns). Present only if
#'   `store_posterior_sample = TRUE`.}
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
#' @details
#' This function is a flexible generalization of \code{HorseTrees}.
#' Instead of using a single Horseshoe prior, it allows specifying different
#' global-local shrinkage configurations for the tree step heights. 
#' Currently, four priors have been implemented. 
#'
#' The \code{horseshoe} prior is the fully Bayesian global-local shrinkage
#' prior, where both the global and local shrinkage parameters are assigned
#' half-Cauchy distributions with scale hyperparameters \code{global_hp} and
#' \code{local_hp}, respectively. The global shrinkage parameter is defined
#' separately for each tree, allowing adaptive regularization per tree.
#'
#' The \code{horseshoe_fw} prior (forest-wide horseshoe) is similar to
#' \code{horseshoe}, except that the global shrinkage parameter is shared
#' across all trees in the forest simultaneously. 
#'
#' The \code{horseshoe_EB} prior is an empirical Bayes variant of the
#' \code{horseshoe} prior. Here, the global shrinkage parameter (\eqn{\tau})
#' is not assigned a prior distribution but instead must be specified directly
#' using \code{global_hp}, while local shrinkage parameters still follow
#' half-Cauchy priors. Note: \eqn{\tau} must be provided by the user; it is 
#' not estimated by the software.
#'
#' The \code{half-cauchy} prior considers only local shrinkage and does not
#' include a global shrinkage component. It places a half-Cauchy prior on each
#' local shrinkage parameter with scale hyperparameter \code{local_hp}.
#' 
#' @examples
#' # Example: Continuous outcome with ShrinkageTrees, two priors
#' n <- 50
#' p <- 3
#' X <- matrix(runif(n * p), ncol = p)
#' X_test <- matrix(runif(n * p), ncol = p)
#' y <- X[, 1] + rnorm(n)
#'
#' # Fit ShrinkageTrees with standard horseshoe prior
#' fit_horseshoe <- ShrinkageTrees(y = y,
#'                                 X_train = X,
#'                                 X_test = X_test,
#'                                 outcome_type = "continuous",
#'                                 number_of_trees = 5,
#'                                 prior_type = "horseshoe",
#'                                 local_hp = 0.1 / sqrt(5),
#'                                 global_hp = 0.1 / sqrt(5),
#'                                 N_post = 10,
#'                                 N_burn = 5,
#'                                 store_posterior_sample = TRUE,
#'                                 verbose = FALSE,
#'                                 seed = 1)
#'
#' # Fit ShrinkageTrees with half-Cauchy prior
#' fit_halfcauchy <- ShrinkageTrees(y = y,
#'                                  X_train = X,
#'                                  X_test = X_test,
#'                                  outcome_type = "continuous",
#'                                  number_of_trees = 5,
#'                                  prior_type = "half-cauchy",
#'                                  local_hp = 1 / sqrt(5),
#'                                  N_post = 10,
#'                                  N_burn = 5,
#'                                  store_posterior_sample = TRUE,
#'                                  verbose = FALSE,
#'                                  seed = 1)
#'
#' # Posterior mean predictions
#' pred_horseshoe <- colMeans(fit_horseshoe$train_predictions_sample)
#' pred_halfcauchy <- colMeans(fit_halfcauchy$train_predictions_sample)
#'
#' # Posteriors of the mean (global average prediction)
#' post_mean_horseshoe <- rowMeans(fit_horseshoe$train_predictions_sample)
#' post_mean_halfcauchy <- rowMeans(fit_halfcauchy$train_predictions_sample)
#'
#' # Posterior mean prediction averages
#' mean_pred_horseshoe <- mean(post_mean_horseshoe)
#' mean_pred_halfcauchy <- mean(post_mean_halfcauchy)
#'
#' @seealso \code{\link{HorseTrees}}, \code{\link{CausalHorseForest}}, \code{\link{CausalShrinkageForest}}
#'
#' @importFrom Rcpp evalCpp
#' @useDynLib ShrinkageTrees, .registration = TRUE
#' @importFrom stats sd qchisq qnorm pnorm runif
#' @export

ShrinkageTrees <- function(y,
                           status = NULL,
                           X_train,
                           X_test = NULL,
                           outcome_type = "continuous",
                           timescale = "time",
                           number_of_trees = 200,
                           prior_type = "horseshoe",
                           local_hp = NULL,
                           global_hp = NULL,
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
                           seed = NULL,
                           verbose = TRUE) {
  
  # Check outcome_type value
  allowed_types <- c("continuous", "binary", "right-censored")
  if (!outcome_type %in% allowed_types) {
    stop("Invalid outcome_type. Please choose 'continuous', 'binary', or 
         'right-censored'.")
  }
  
  # Check prior_type value
  allowed_prior <- c("horseshoe", "horseshoe_fw", "horseshoe_EB", "half-cauchy")
  if (!prior_type %in% allowed_prior) {
    stop("Invalid prior_type. Choose 'horseshoe', 'horseshoe_fw', 'horseshoe_EB', or 'half-cauchy'.")
  }
  
  # Prior-specific checks
  if (prior_type %in% c("horseshoe", "horseshoe_fw", "horseshoe_EB")) {
    if (is.null(local_hp) || is.null(global_hp)) {
      stop("For prior_type = 'horseshoe', 'horseshoe_fw', or 'horseshoe_EB', you must provide both local_hp and global_hp.")
    }
  }
  
  if (prior_type == "half-cauchy") {
    if (is.null(local_hp)) {
      stop("For prior_type = 'half-cauchy', you must provide local_hp.")
    }
    if (!is.null(global_hp)) {
      warning("global_hp is ignored for 'half-cauchy'. If you want to fix the global parameter, consider using 'horseshoe_EB'.")
    }
    global_hp <- 1
    prior_type <- "halfcauchy"
  }
  
  if (prior_type == "horseshoe_EB") prior_type <- "halfcauchy"

  # Check consistency with status argument
  if (outcome_type == "right-censored" && is.null(status)) {
    stop("You specified outcome_type = 'right-censored', but did not provide a 
         'status' vector.")
  }
  
  if (outcome_type != "right-censored" && !is.null(status)) {
    warning("You provided a 'status' vector, but outcome_type is not 
            'right-censored'. The 'status' vector will be ignored.")
  }
  
  # Check binary data consistency
  if (outcome_type != "binary" && all(y %in% c(0, 1))) {
    warning("The outcome y contains only 0 and 1, but outcome_type is not set to
            'binary'. Consider setting outcome_type = 'binary'.")
  }
  
  # Check survival data and timescale
  if (outcome_type == "right-censored" && timescale == "time" && any(y < 0)) {
    stop("Outcome contains negative values, but timescale = 'time' for survival 
         data requires non-negative times.")
  }
  
  # Retrieve dimensions of training data
  n_train <- nrow(X_train)
  p_features <- ncol(X_train)
  
  # Check if dimensions match with test data
  if (!is.null(X_test)) {
    n_test <- nrow(X_test)
    if (ncol(X_test) != p_features) {
      stop("The number of covariates in X_test must match X_train.")
    }
    X_test <- as.numeric(t(X_test))
  } else {
    n_test <- 1
    X_test <- as.numeric(colMeans(X_train))
  }
  
  # Force data types to be numeric and plain arrays
  N_post <- as.integer(N_post)[1]
  N_burn <- as.integer(N_burn)[1]
  power <- as.numeric(power)[1]
  base <- as.numeric(base)[1]
  p_grow <- as.numeric(p_grow)[1]
  p_prune <- as.numeric(p_prune)[1]
  X_train <- as.numeric(t(X_train))
  
  # Set a random seed if not provided
  # By taking a random number, we ensure compatibility with set.seed()
  if (is.null(seed)) seed <- as.integer(runif(1, 1, 1000000))
  
  if (outcome_type == "right-censored") {
    
    # Convert y to numeric for C++ compatibility
    y <- as.numeric(y)
    
    # Log-transform survival times if timescale = "time"
    if (timescale == "time") {
      y <- log(y)
    }
    
    # Obtain estimated mean and standard deviation using censored_info()
    cens_inf <- censored_info(y, status)
    
    # Center the data using estimated mean
    y_mean <- cens_inf$mu
    y <- y - y_mean
    
    # Determine sigma (timescale parameter) and whether it is known
    if (is.null(sigma)) {
      sigma_hat <- cens_inf$sd
      sigma_known <- FALSE
    } else {
      sigma_hat <- sigma
      sigma_known <- TRUE
    }
    
    # Standardize the centered data
    y <- y / sigma_hat
    
    # Set survival flag
    survival <- TRUE
    
    # Compute lambda parameter for error distribution prior
    qchi <- qchisq(1.0 - q, nu)
    lambda <- (sigma_hat^2 * qchi) / nu
    
    # Fit a HorseTrees model
    fit <- HorseTrees_cpp(
      nSEXP = n_train,
      pSEXP = p_features,
      n_testSEXP = n_test,
      X_trainSEXP = X_train,
      ySEXP = y,
      status_indicatorSEXP = status,
      is_survivalSEXP = survival,
      X_testSEXP = X_test,
      number_of_treesSEXP = number_of_trees,
      N_postSEXP = N_post,
      N_burnSEXP = N_burn,
      delayed_proposalSEXP = delayed_proposal,
      powerSEXP = power,
      baseSEXP = base,
      p_growSEXP = p_grow,
      p_pruneSEXP = p_prune,
      nuSEXP = nu,
      lambdaSEXP = lambda,
      sigmaSEXP = sigma_hat,
      sigma_knownSEXP = sigma_known,
      omegaSEXP = 1,
      param1SEXP = local_hp,
      param2SEXP = global_hp,
      prior_typeSEXP = prior_type,
      reversibleSEXP = TRUE,
      store_parametersSEXP = FALSE,
      store_posterior_sampleSEXP = store_posterior_sample,
      n1SEXP = seed,
      n2SEXP = 420,
      verboseSEXP = verbose
    )
    
    if (timescale == "time") {
      fit$train_predictions <- exp(fit$train_predictions * sigma_hat + y_mean)
      fit$test_predictions <- exp(fit$test_predictions * sigma_hat + y_mean)
      if (store_posterior_sample) {
        fit$train_predictions_sample <- exp(fit$train_predictions_sample * sigma_hat + y_mean)
        fit$test_predictions_sample <- exp(fit$test_predictions_sample * sigma_hat + y_mean)
      }
    } else { # If timescale is "log"
      fit$train_predictions <- fit$train_predictions * sigma_hat + y_mean
      fit$test_predictions <- fit$test_predictions * sigma_hat + y_mean
      if (store_posterior_sample) { 
        fit$train_predictions_sample <- fit$train_predictions_sample * sigma_hat + y_mean
        fit$test_predictions_sample <- fit$test_predictions_sample * sigma_hat + y_mean
      }
    }
    
    # If binary
  } else if (outcome_type == "binary") {
    y <- as.numeric(y)
    latent_threshold <- qnorm(mean(y))
    
    fit <- probitHorseTrees_cpp(
      nSEXP = n_train,
      pSEXP = p_features,
      n_testSEXP = n_test,
      X_trainSEXP = X_train,
      ySEXP = y,
      X_testSEXP = X_test,
      number_of_treesSEXP = number_of_trees,
      N_postSEXP = N_post,
      N_burnSEXP = N_burn,
      delayed_proposalSEXP = delayed_proposal,
      powerSEXP = power,
      baseSEXP = base,
      p_growSEXP = p_grow,
      p_pruneSEXP = p_prune,
      omegaSEXP = 1,
      latent_thresholdSEXP = latent_threshold,
      param1SEXP = local_hp,
      param2SEXP = global_hp,
      prior_typeSEXP = prior_type,
      reversibleSEXP = TRUE,
      store_posterior_sampleSEXP = store_posterior_sample,
      n1SEXP = seed,
      n2SEXP = 420,
      verboseSEXP = verbose
    )
    
    # Add the pnorm transform
    fit$train_probabilities <- pnorm(fit$train_predictions)
    fit$test_probabilities<- pnorm(fit$test_predictions)
    if (store_posterior_sample) {
      fit$train_probabilities_sample <- pnorm(fit$train_predictions_sample)
      fit$test_probabilities_sample <- pnorm(fit$test_predictions_sample)
    }
    
    # Otherwise, continuous
  } else {
    
    # Force outcome to plain numeric vector
    y <- as.numeric(y)
    
    # Create dummy status vector (not used for continuous)
    status <- rep(1, n_train)
    survival <- FALSE
    
    # Determine prior guess of sigma 
    if (is.null(sigma)) {
      sigma_hat <- sd(y)      # Estimate sigma from data
      sigma_known <- FALSE
    } else {
      sigma_hat <- sigma      # Use provided sigma
      sigma_known <- TRUE
    }
    
    # Compute hyperparameters of error variance prior
    qchi <- qchisq(1.0 - q, nu)
    lambda <- (sigma_hat^2 * qchi) / nu
    
    # Compute mean and standardize  the outcome
    y_mean <- mean(y)
    y <- y - y_mean
    y <- y / sigma_hat  # Standardize y
    
    fit <- HorseTrees_cpp(
      nSEXP = n_train,
      pSEXP = p_features,
      n_testSEXP = n_test,
      X_trainSEXP = X_train,
      ySEXP = y,
      status_indicatorSEXP = status,
      is_survivalSEXP = survival,
      X_testSEXP = X_test,
      number_of_treesSEXP = number_of_trees,
      N_postSEXP = N_post,
      N_burnSEXP = N_burn,
      delayed_proposalSEXP = delayed_proposal,
      powerSEXP = power,
      baseSEXP = base,
      p_growSEXP = p_grow,
      p_pruneSEXP = p_prune,
      nuSEXP = nu,
      lambdaSEXP = lambda,
      sigmaSEXP = sigma_hat,
      sigma_knownSEXP = sigma_known,
      omegaSEXP = 1,
      param1SEXP = local_hp,
      param2SEXP = global_hp,
      prior_typeSEXP = prior_type,
      reversibleSEXP = TRUE,
      store_parametersSEXP = FALSE,
      store_posterior_sampleSEXP = store_posterior_sample,
      n1SEXP = seed,
      n2SEXP = 420,
      verboseSEXP = verbose
    )
    
    fit$train_predictions <- fit$train_predictions * sigma_hat + y_mean
    fit$test_predictions <- fit$test_predictions * sigma_hat + y_mean
    if (store_posterior_sample) {
      fit$train_predictions_sample <- fit$train_predictions_sample * sigma_hat + y_mean
      fit$test_predictions_sample <- fit$test_predictions_sample * sigma_hat + y_mean
    }
  }
  
  return(fit)
}
