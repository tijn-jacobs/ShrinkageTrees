#' Horseshoe Prior Bayesian Tree Model (HorseTrees)
#'
#' Fits a Bayesian Horseshoe Trees model. 
#' Supports continuous, binary, and right-censored (survival) outcomes.
#'
#' For continuous outcomes, the model centers and optionally standardizes the outcome using a prior guess of the standard deviation.
#' For binary outcomes, the function uses a probit link formulation.
#' For right-censored outcomes (survival data), the function can handle follow-up times either on the original time scale or log-transformed.
#'
#' @param y Outcome vector. Numeric. Can represent continuous outcomes, binary outcomes (0/1), or follow-up times for survival data.
#' @param X_train Covariate matrix for training. Each row corresponds to an observation, and each column to a covariate.
#' @param X_test Optional covariate matrix for test data. If NULL, defaults to the mean of the training covariates.
#' @param status Optional censoring indicator vector (1 = event occurred, 0 = censored). Required if `outcome_type = "right-censored"`.
#' @param outcome_type Type of outcome. One of `"continuous"`, `"binary"`, or `"right-censored"`.
#' @param timescale Indicates the scale of follow-up times. Options are `"time"` (nonnegative follow-up times, will be log-transformed internally) 
#'   or `"log"` (already log-transformed). Only used when `outcome_type = "right-censored"`.
#' @param number_of_trees Number of trees in the ensemble. Default is 200.
#' @param k Horseshoe hyperparameter (default 0.1).
#' @param power Power parameter for tree structure prior. Default is 2.0.
#' @param base Base parameter for tree structure prior. Default is 0.95.
#' @param p_grow Probability of proposing a grow move. Default is 0.4.
#' @param p_prune Probability of proposing a prune move. Default is 0.4.
#' @param nu Degrees of freedom for the error distribution prior. Default is 3.
#' @param q Quantile hyperparameter for the error variance prior. Default is 0.90.
#' @param sigma Optional known value for error standard deviation. If NULL, estimated from data.
#' @param N_post Number of posterior samples to store. Default is 1000.
#' @param N_burn Number of burn-in iterations. Default is 1000.
#' @param delayed_proposal Number of delayed iterations before proposal. Only for reversible updates. Default is 5.
#' @param store_posterior_sample Logical; whether to store posterior samples for each iteration. Default is TRUE.
#' @param seed Random seed for reproducibility. Default uses system time.
#' @param verbose Logical; whether to print verbose output. Default is TRUE.
#'
#' @return A list containing fitted object with posterior mean predictions, test predictions (if provided), and optionally stored posterior samples.
#' @examples
#' # Example usage with continuous outcome
#' set.seed(1)
#' X <- matrix(rnorm(100 * 5), ncol = 5)
#' y <- X[,1] + rnorm(100)
#' fit <- HorseTrees(y, X, outcome_type = "continuous")
#' 
#' # Example usage with binary outcome
#' y_bin <- ifelse(X[,1] + rnorm(100) > 0, 1, 0)
#' fit_bin <- HorseTrees(y_bin, X, outcome_type = "binary", verbose = FALSE)
#' 
#' # Example usage with right-censored outcome (survival data)
#' library(survival)
#' time <- rexp(100, rate = 0.1)
#' status <- rbinom(100, 1, 0.7)
#' fit_surv <- HorseTrees(time, X, status = status, 
#'                        outcome_type = "right-censored", timescale = "time", 
#'                        verbose = FALSE)
#'
#' @useDynLib ShrinkageTrees, .registration = TRUE
#' @importFrom stats sd qchisq qnorm
#' @export

HorseTrees <- function(y,
                       X_train,
                       X_test = NULL,
                       status = NULL,
                       outcome_type = "continuous",
                       timescale = "time",
                       number_of_trees = 200,
                       k = 0.1,
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
    stop("Invalid outcome_type. Please choose 'continuous', 'binary', or 'right-censored'.")
  }

  # Check consistency with status argument
  if (outcome_type == "right-censored" && is.null(status)) {
    stop("You specified outcome_type = 'right-censored', but did not provide a 'status' vector.")
  }

  if (outcome_type != "right-censored" && !is.null(status)) {
    warning("You provided a 'status' vector, but outcome_type is not 'right-censored'. The 'status' vector will be ignored.")
  }

  # Check binary data consistency
  if (outcome_type != "binary" && all(y %in% c(0, 1))) {
    warning("The outcome y contains only 0 and 1, but outcome_type is not set to 'binary'. Consider setting outcome_type = 'binary'.")
  }

  # Check survival data and timescale
  if (outcome_type == "right-censored" && timescale == "time" && any(y < 0)) {
    stop("Outcome contains negative values, but timescale = 'time' for survival data requires non-negative times.")
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
  if (is.null(seed)) seed <- as.integer(Sys.time())

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
      param1SEXP = k / sqrt(number_of_trees),
      param2SEXP = k / sqrt(number_of_trees),
      prior_typeSEXP = "horseshoe",
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
      fit$train_predictions_sample <- exp(fit$train_predictions_sample * sigma_hat + y_mean)
      fit$test_predictions_sample <- exp(fit$test_predictions_sample * sigma_hat + y_mean)
    } else { # If timescale is "log"
      fit$train_predictions <- fit$train_predictions * sigma_hat + y_mean
      fit$test_predictions <- fit$test_predictions * sigma_hat + y_mean
      fit$train_predictions_sample <- fit$train_predictions_sample * sigma_hat + y_mean
      fit$test_predictions_sample <- fit$test_predictions_sample * sigma_hat + y_mean
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
      param1SEXP = k / sqrt(number_of_trees),
      param2SEXP = k / sqrt(number_of_trees),
      prior_typeSEXP = "horseshoe",
      reversibleSEXP = TRUE,
      store_posterior_sampleSEXP = store_posterior_sample,
      n1SEXP = seed,
      n2SEXP = 420,
      verboseSEXP = verbose
    )

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
      param1SEXP = k / sqrt(number_of_trees),
      param2SEXP = k / sqrt(number_of_trees),
      prior_typeSEXP = "horseshoe",
      reversibleSEXP = TRUE,
      store_parametersSEXP = FALSE,
      store_posterior_sampleSEXP = store_posterior_sample,
      n1SEXP = seed,
      n2SEXP = 420,
      verboseSEXP = verbose
    )
    
    fit$train_predictions <- fit$train_predictions * sigma_hat + y_mean
    fit$test_predictions <- fit$test_predictions * sigma_hat + y_mean
    fit$train_predictions_sample <- fit$train_predictions_sample * sigma_hat + y_mean
    fit$test_predictions_sample <- fit$test_predictions_sample * sigma_hat + y_mean
  }
  
  return(fit)
}
