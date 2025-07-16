#' General Causal Shrinkage Forests
#'
#' Fits a (Bayesian) Causal Shrinkage Forest model for estimating heterogeneous treatment effects.
#' This function generalizes \code{\link{CausalHorseForest}} by allowing flexible
#' global-local shrinkage priors on the step heights in both the control and treatment forests.
#' It supports continuous and right-censored survival outcomes.
#'
#' @param y Outcome vector. Numeric. Represents continuous outcomes or follow-up times.
#' @param status Optional event indicator vector (1 = event occurred, 0 = censored). 
#' Required when \code{outcome_type = "right-censored"}.
#' @param X_train_control Covariate matrix for the control forest. Rows correspond to samples,
#' columns to covariates.
#' @param X_train_treat Covariate matrix for the treatment forest.
#' @param treatment_indicator_train Vector indicating treatment assignment for training samples 
#' (1 = treated, 0 = control).
#' @param X_test_control Optional covariate matrix for control forest test data. Defaults to 
#' column means of \code{X_train_control} if NULL.
#' @param X_test_treat Optional covariate matrix for treatment forest test data. Defaults to 
#' column means of \code{X_train_treat} if NULL.
#' @param treatment_indicator_test Optional vector indicating treatment assignment for test data.
#' @param outcome_type Type of outcome: one of \code{"continuous"} or \code{"right-censored"}. 
#' Default is \code{"continuous"}.
#' @param timescale For survival outcomes: either \code{"time"} (original scale, log-transformed 
#' internally) or \code{"log"} (already log-transformed). Default is \code{"time"}.
#' @param number_of_trees_control Number of trees in the control forest. Default is 200.
#' @param number_of_trees_treat Number of trees in the treatment forest. Default is 200.
#' @param prior_type_control Type of prior on control forest step heights. One of 
#' \code{"horseshoe"}, \code{"horseshoe_fw"}, \code{"horseshoe_EB"}, or \code{"half-cauchy"}. 
#' Default is \code{"horseshoe"}.
#' @param prior_type_treat Type of prior on treatment forest step heights. Same options as 
#' \code{prior_type_control}.
#' @param local_hp_control Local hyperparameter controlling shrinkage on individual steps 
#' (control forest). Required for all prior types.
#' @param local_hp_treat Local hyperparameter for treatment forest.
#' @param global_hp_control Global hyperparameter for control forest. Required for horseshoe-type
#' priors; ignored for \code{"half-cauchy"}.
#' @param global_hp_treat Global hyperparameter for treatment forest.
#' @param power Power parameter for tree structure prior. Default is 2.0.
#' @param base Base parameter for tree structure prior. Default is 0.95.
#' @param p_grow Probability of proposing a grow move. Default is 0.4.
#' @param p_prune Probability of proposing a prune move. Default is 0.4.
#' @param nu Degrees of freedom for the error variance prior. Default is 3.
#' @param q Quantile parameter for error variance prior. Default is 0.90.
#' @param sigma Optional known standard deviation of the outcome. If NULL, estimated from data.
#' @param N_post Number of posterior samples to store. Default is 5000.
#' @param N_burn Number of burn-in iterations. Default is 5000.
#' @param delayed_proposal Number of delayed iterations before proposal updates. Default is 5.
#' @param store_posterior_sample Logical; whether to store posterior samples of predictions. 
#' Default is \code{FALSE}.
#' @param seed Random seed for reproducibility. Default is \code{NULL}.
#' @param verbose Logical; whether to print verbose output. Default is \code{TRUE}.
#'
#' @return A list containing:
#' \describe{
#'   \item{train_predictions}{Posterior mean predictions on training data (combined forest).}
#'   \item{test_predictions}{Posterior mean predictions on test data (combined forest).}
#'   \item{train_predictions_control}{Estimated control outcomes on training data.}
#'   \item{test_predictions_control}{Estimated control outcomes on test data.}
#'   \item{train_predictions_treat}{Estimated treatment effects on training data.}
#'   \item{test_predictions_treat}{Estimated treatment effects on test data.}
#'   \item{sigma}{Vector of posterior samples for the error standard deviation.}
#'   \item{acceptance_ratio_control}{Average acceptance ratio in control forest.}
#'   \item{acceptance_ratio_treat}{Average acceptance ratio in treatment forest.}
#'   \item{train_predictions_sample_control}{Matrix of posterior samples for control predictions 
#'   (if \code{store_posterior_sample = TRUE}).}
#'   \item{test_predictions_sample_control}{Matrix of posterior samples for control predictions 
#'   (if \code{store_posterior_sample = TRUE}).}
#'   \item{train_predictions_sample_treat}{Matrix of posterior samples for treatment effects 
#'   (if \code{store_posterior_sample = TRUE}).}
#'   \item{test_predictions_sample_treat}{Matrix of posterior samples for treatment effects 
#'   (if \code{store_posterior_sample = TRUE}).}
#' }
#'
#' @details
#' This function is a flexible generalization of \code{CausalHorseForest}.
#' The Causal Shrinkage Forest model decomposes the outcome into a prognostic 
#' (control) and a treatment effect part. Each part is modeled by its own 
#' shrinkage tree ensemble, with separate flexible global-local shrinkage 
#' priors. It is particularly useful for estimating heterogeneous treatment 
#' effects in high-dimensional settings.
#'
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
#' # Example: Continuous outcome, homogenuous treatment effect, two priors
#' n <- 50
#' p <- 3
#' X <- matrix(runif(n * p), ncol = p)
#' X_treat <- X_control <- X
#' treat <- rbinom(n, 1, X[,1])
#' tau <- 2
#' y <- X[, 1] + (0.5 - treat) * tau + rnorm(n)
#' 
#' # Fit a standard Causal Horseshoe Forest
#' fit_horseshoe <- CausalShrinkageForest(y = y,
#'                                        X_train_control = X_control,
#'                                        X_train_treat = X_treat,
#'                                        treatment_indicator_train = treat,
#'                                        outcome_type = "continuous",
#'                                        number_of_trees_treat = 5,
#'                                        number_of_trees_control = 5,
#'                                        prior_type_control = "horseshoe",
#'                                        prior_type_treat = "horseshoe",
#'                                        local_hp_control = 0.1/sqrt(5),
#'                                        local_hp_treat = 0.1/sqrt(5),
#'                                        global_hp_control = 0.1/sqrt(5),
#'                                        global_hp_treat = 0.1/sqrt(5),
#'                                        N_post = 10,
#'                                        N_burn = 5,
#'                                        store_posterior_sample = TRUE,
#'                                        verbose = FALSE,
#'                                        seed = 1
#' )
#' 
#' # Fit a Causal Shrinkage Forest with half-cauchy prior
#' fit_halfcauchy <- CausalShrinkageForest(y = y,
#'                                         X_train_control = X_control,
#'                                         X_train_treat = X_treat,
#'                                         treatment_indicator_train = treat,
#'                                         outcome_type = "continuous",
#'                                         number_of_trees_treat = 5,
#'                                         number_of_trees_control = 5,
#'                                         prior_type_control = "half-cauchy",
#'                                         prior_type_treat = "half-cauchy",
#'                                         local_hp_control = 1/sqrt(5),
#'                                         local_hp_treat = 1/sqrt(5),
#'                                         N_post = 10,
#'                                         N_burn = 5,
#'                                         store_posterior_sample = TRUE,
#'                                         verbose = FALSE,
#'                                         seed = 1
#' )
#' 
#' # Posterior mean CATEs
#' CATE_horseshoe <- colMeans(fit_horseshoe$train_predictions_sample_treat)
#' CATE_halfcauchy <- colMeans(fit_halfcauchy$train_predictions_sample_treat)
#' 
#' # Posteriors of the ATE
#' post_ATE_horseshoe <- rowMeans(fit_horseshoe$train_predictions_sample_treat)
#' post_ATE_halfcauchy <- rowMeans(fit_halfcauchy$train_predictions_sample_treat)
#' 
#' # Posterior mean ATE
#' ATE_horseshoe <- mean(post_ATE_horseshoe)
#' ATE_halfcauchy <- mean(post_ATE_halfcauchy)
#' 
#' 
#' @seealso \code{\link{CausalHorseForest}}, \code{\link{ShrinkageTrees}}, 
#' \code{\link{HorseTrees}}
#'
#' @importFrom Rcpp evalCpp
#' @useDynLib ShrinkageTrees, .registration = TRUE
#' @importFrom stats sd qchisq qnorm runif coef
#' @export
CausalShrinkageForest <- function(y,
                                  status = NULL,
                                  X_train_control,
                                  X_train_treat,
                                  treatment_indicator_train,
                                  X_test_control = NULL,
                                  X_test_treat = NULL,
                                  treatment_indicator_test = NULL,
                                  outcome_type = "continuous",
                                  timescale = "time",
                                  number_of_trees_control = 200,
                                  number_of_trees_treat = 200,
                                  prior_type_control = "horseshoe",
                                  prior_type_treat = "horseshoe",
                                  local_hp_control,
                                  local_hp_treat,
                                  global_hp_control = NULL,
                                  global_hp_treat = NULL,
                                  power = 2.0,
                                  base = 0.95,
                                  p_grow = 0.4,
                                  p_prune = 0.4,
                                  nu = 3,
                                  q = 0.90,
                                  sigma = NULL,
                                  N_post = 5000,
                                  N_burn = 5000,
                                  delayed_proposal = 5,
                                  store_posterior_sample = FALSE,
                                  seed = NULL,
                                  verbose = TRUE) {

  
  # Check prior_type value
  allowed_prior <- c("horseshoe", "horseshoe_fw", "horseshoe_EB", "half-cauchy")
  if (!prior_type_control %in% allowed_prior) {
    stop("Invalid prior_type_control Choose 'horseshoe', 'horseshoe_fw', 
         'horseshoe_EB', or 'half-cauchy'.")
  }
  
  if (!prior_type_treat %in% allowed_prior) {
    stop("Invalid prior_type_treat. Choose 'horseshoe', 'horseshoe_fw', 
         'horseshoe_EB', or 'half-cauchy'.")
  }
  
  # Prior-specific checks
  if (prior_type_control %in% c("horseshoe", "horseshoe_fw", "horseshoe_EB")) {
    if (is.null(local_hp_control) || is.null(global_hp_control)) {
      stop("For prior_type_control = 'horseshoe', 'horseshoe_fw', or 
           'horseshoe_EB', you must provide both local_hp and global_hp.")
    }
  }
  
  if (prior_type_treat %in% c("horseshoe", "horseshoe_fw", "horseshoe_EB")) {
    if (is.null(local_hp_treat) || is.null(global_hp_treat)) {
      stop("For prior_type_treat = 'horseshoe', 'horseshoe_fw', or 
           'horseshoe_EB', you must provide both local_hp and global_hp.")
    }
  }
  
  if (prior_type_control == "half-cauchy") {
    if (is.null(local_hp_control)) {
      stop("For prior_type = 'half-cauchy', you must provide local_hp_control.")
    }
    if (!is.null(global_hp_control)) {
      warning("global_hp_control is ignored for 'half-cauchy'. If you want to 
              fix the global parameter, consider using 'horseshoe_EB'.")
    }
    global_hp_control <- 1
    prior_type_control <- "halfcauchy"
  }
  
  if (prior_type_treat == "half-cauchy") {
    if (is.null(local_hp_treat)) {
      stop("For prior_type_treat = 'half-cauchy', you must provide 
           local_hp_treat.")
    }
    if (!is.null(global_hp_treat)) {
      warning("global_hp_treat is ignored for 'half-cauchy'. If you want to fix 
              the global parameter, consider using 'horseshoe_EB'.")
    }
    global_hp_treat <- 1
    prior_type_treat <- "halfcauchy"
  }
  
  if (prior_type_control == "horseshoe_EB") prior_type_control <- "halfcauchy"
  if (prior_type_treat == "horseshoe_EB") prior_type_treat <- "halfcauchy"
  
  # Check outcome_type value
  allowed_types <- c("continuous", "right-censored")
  if (!outcome_type %in% allowed_types) {
    stop("Invalid outcome_type. Please choose 'continuous' or 
         'right-censored'.")
  }
  
  # Check consistency with status argument
  if (outcome_type == "right-censored" && is.null(status)) {
    stop("You specified outcome_type = 'right-censored', but did not provide a 
         'status' vector.")
  }
  
  if (outcome_type != "right-censored" && !is.null(status)) {
    warning("You provided a 'status' vector, but outcome_type is not 
            'right-censored'. The 'status' vector will be ignored.")
  }
  
  # Check survival data and timescale
  if (outcome_type == "right-censored" && timescale == "time" && any(y < 0)) {
    stop("Outcome contains negative values, but timescale = 'time' for survival 
         data requires non-negative times.")
  }
  
  # Retrieve dimensions of training data
  n_train <- nrow(X_train_control)
  p_control <- ncol(X_train_control)
  p_treat   <- ncol(X_train_treat)
  
  # Check matching row numbers
  if (nrow(X_train_control) != length(y)) {
    stop("X_train_control rows must match length of y.")
  }
  if (nrow(X_train_treat) != length(y)) {
    stop("X_train_treat rows must match length of y.")
  }
  if (length(treatment_indicator_train) != length(y)) {
    stop("treatment_indicator_train must match length of y.")
  }
  
  # If test provided, check
  if (!is.null(X_test_control) && !is.null(X_test_treat)) {
    n_test <- nrow(X_test_control)
    if (ncol(X_test_control) != p_control || ncol(X_test_treat) != p_treat) {
      stop("Number of columns in X_test_control or X_test_treat does not match 
           training data.")
    }
    if (!is.null(treatment_indicator_test) && 
        length(treatment_indicator_test) != n_test) {
      stop("treatment_indicator_test length must match number of test rows.")
    }
    X_test_control <- as.numeric(t(X_test_control))
    X_test_treat <- as.numeric(t(X_test_treat))
  } else {
    n_test <- 1
    X_test_control <- as.numeric(colMeans(X_train_control))
    X_test_treat <- as.numeric(colMeans(X_train_treat))
    treatment_indicator_test <- as.integer(rep(1, n_test))
  }
  
  # Treatment indicators as integer
  treatment_indicator_train <- as.integer(treatment_indicator_train)
  if (is.null(treatment_indicator_test)) {
    treatment_indicator_test <- as.integer(rep(1, n_test))
  } else {
    treatment_indicator_test <- as.integer(treatment_indicator_test)
  }
  
  # Force data types to be numeric and plain arrays
  N_post <- as.integer(N_post)[1]
  N_burn <- as.integer(N_burn)[1]
  power <- as.numeric(power)[1]
  base <- as.numeric(base)[1]
  p_grow <- as.numeric(p_grow)[1]
  p_prune <- as.numeric(p_prune)[1]
  X_train_treat <- as.numeric(t(X_train_treat))
  X_train_control <- as.numeric(t(X_train_control))
  
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
    
    # Fit a Causal Horseshoe Forest
    fit <- CausalHorseForest_cpp(
      nSEXP = n_train,
      p_treatSEXP = p_treat,
      p_controlSEXP = p_control,
      X_train_treatSEXP = X_train_treat,
      X_train_controlSEXP = X_train_control,
      ySEXP = y,
      status_indicatorSEXP = status,
      is_survivalSEXP = survival,
      treatment_indicatorSEXP = treatment_indicator_train,
      n_testSEXP = n_test,
      X_test_controlSEXP = X_test_control,
      X_test_treatSEXP = X_test_treat,
      treatment_indicator_testSEXP = treatment_indicator_test,
      no_trees_treatSEXP = number_of_trees_treat,
      power_treatSEXP = power,
      base_treatSEXP = base,
      p_grow_treatSEXP = p_grow,
      p_prune_treatSEXP = p_prune,
      omega_treatSEXP = 1/2,
      prior_type_treatSEXP = prior_type_treat,
      param1_treatSEXP = local_hp_treat,
      param2_treatSEXP = global_hp_treat,
      reversible_treatSEXP = TRUE,
      no_trees_controlSEXP = number_of_trees_control,
      power_controlSEXP = power,
      base_controlSEXP = base,
      p_grow_controlSEXP = p_grow,
      p_prune_controlSEXP = p_prune,
      omega_controlSEXP = 1/2,
      prior_type_controlSEXP = prior_type_control,
      param1_controlSEXP = local_hp_control,
      param2_controlSEXP = global_hp_control,
      reversible_controlSEXP = TRUE,
      sigma_knownSEXP = sigma_known,
      sigmaSEXP = sigma_hat,
      lambdaSEXP = lambda,
      nuSEXP = nu,
      N_postSEXP = N_post,
      N_burnSEXP = N_burn,
      delayed_proposalSEXP = delayed_proposal,
      store_parametersSEXP = FALSE,
      max_stored_leafsSEXP = 1,
      store_posterior_sample_controlSEXP = store_posterior_sample,
      store_posterior_sample_treatSEXP = store_posterior_sample,
      n1SEXP = seed,
      n2SEXP = 420,
      verboseSEXP = verbose
    )
    
    if (timescale == "time") {
      
      # Total
      fit$train_predictions <- exp(fit$train_predictions * sigma_hat + y_mean)
      fit$test_predictions <- exp(fit$test_predictions * sigma_hat + y_mean)
      
      # Control
      fit$train_predictions_control <- exp(fit$train_predictions_control * 
                                             sigma_hat + y_mean)
      fit$test_predictions_control <- exp(fit$test_predictions_control * 
                                            sigma_hat + y_mean)
      
      # Treatment
      fit$train_predictions_treat <- exp(fit$train_predictions_treat * 
                                           sigma_hat)
      fit$test_predictions_treat <- exp(fit$test_predictions_treat * sigma_hat)
      
      # Posterior samples
      if (store_posterior_sample) { #
        fit$train_predictions_sample_control <- 
          exp(fit$train_predictions_sample_control * sigma_hat + y_mean)
        fit$test_predictions_sample_control <- 
          exp(fit$test_predictions_sample_control * sigma_hat + y_mean)
        fit$train_predictions_sample_treat <- 
          exp(fit$train_predictions_sample_treat * sigma_hat)
        fit$test_predictions_sample_treat <- 
          exp(fit$test_predictions_sample_treat * sigma_hat)
      }
    } else {
      # Total
      fit$train_predictions <- fit$train_predictions * sigma_hat + y_mean
      fit$test_predictions <- fit$test_predictions * sigma_hat + y_mean
      
      # Control
      fit$train_predictions_control <- 
        fit$train_predictions_control * sigma_hat + y_mean
      fit$test_predictions_control <- 
        fit$test_predictions_control * sigma_hat + y_mean
      
      # Treatment
      fit$train_predictions_treat <- fit$train_predictions_treat * sigma_hat
      fit$test_predictions_treat <- fit$test_predictions_treat * sigma_hat
      
      # Posterior samples
      if (store_posterior_sample) {
        fit$train_predictions_sample_control <- 
          fit$train_predictions_sample_control * sigma_hat + y_mean
        fit$test_predictions_sample_control <- 
          fit$test_predictions_sample_control * sigma_hat + y_mean
        fit$train_predictions_sample_treat <- 
          fit$train_predictions_sample_treat * sigma_hat
        fit$test_predictions_sample_treat <- 
          fit$test_predictions_sample_treat * sigma_hat
      }
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
    
    # Fit a Causal Horseshoe Forest
    fit <- CausalHorseForest_cpp(
      nSEXP = n_train,
      p_treatSEXP = p_treat,
      p_controlSEXP = p_control,
      X_train_treatSEXP = X_train_treat,
      X_train_controlSEXP = X_train_control,
      ySEXP = y,
      status_indicatorSEXP = status,
      is_survivalSEXP = survival,
      treatment_indicatorSEXP = treatment_indicator_train,
      n_testSEXP = n_test,
      X_test_controlSEXP = X_test_control,
      X_test_treatSEXP = X_test_treat,
      treatment_indicator_testSEXP = treatment_indicator_test,
      no_trees_treatSEXP = number_of_trees_treat,
      power_treatSEXP = power,
      base_treatSEXP = base,
      p_grow_treatSEXP = p_grow,
      p_prune_treatSEXP = p_prune,
      omega_treatSEXP = 1/2,
      prior_type_treatSEXP = prior_type_treat,
      param1_treatSEXP = local_hp_treat,
      param2_treatSEXP = global_hp_treat,
      reversible_treatSEXP = TRUE,
      no_trees_controlSEXP = number_of_trees_control,
      power_controlSEXP = power,
      base_controlSEXP = base,
      p_grow_controlSEXP = p_grow,
      p_prune_controlSEXP = p_prune,
      omega_controlSEXP = 1/2,
      prior_type_controlSEXP = prior_type_control,
      param1_controlSEXP = local_hp_control,
      param2_controlSEXP = global_hp_control,
      reversible_controlSEXP = TRUE,
      sigma_knownSEXP = sigma_known,
      sigmaSEXP = sigma_hat,
      lambdaSEXP = lambda,
      nuSEXP = nu,
      N_postSEXP = N_post,
      N_burnSEXP = N_burn,
      delayed_proposalSEXP = delayed_proposal,
      store_parametersSEXP = FALSE,
      max_stored_leafsSEXP = 1,
      store_posterior_sample_controlSEXP = store_posterior_sample,
      store_posterior_sample_treatSEXP = store_posterior_sample,
      n1SEXP = seed,
      n2SEXP = 420,
      verboseSEXP = verbose
    )
    
    # Total
    fit$train_predictions <- fit$train_predictions * sigma_hat + y_mean
    fit$test_predictions <- fit$test_predictions * sigma_hat + y_mean
    
    # Control
    fit$train_predictions_control <- fit$train_predictions_control * 
      sigma_hat + y_mean
    fit$test_predictions_control <- fit$test_predictions_control * 
      sigma_hat + y_mean
    
    # Treatment
    fit$train_predictions_treat <- fit$train_predictions_treat * sigma_hat 
    fit$test_predictions_treat <- fit$test_predictions_treat * sigma_hat 
    
    # Posterior samples
    if (!is.null(fit$train_predictions_sample)) {
      fit$train_predictions_sample <- 
        fit$train_predictions_sample *  sigma_hat + y_mean
      fit$test_predictions_sample <- 
        fit$test_predictions_sample * sigma_hat +  y_mean
      fit$train_predictions_sample_control <- 
        fit$train_predictions_sample_control * sigma_hat + y_mean
      fit$test_predictions_sample_control <- 
        fit$test_predictions_sample_control * sigma_hat + y_mean
      fit$train_predictions_sample_treat <- 
        fit$train_predictions_sample_treat * sigma_hat 
      fit$test_predictions_sample_treat <- 
        fit$test_predictions_sample_treat * sigma_hat 
    }
    
  }
  
  return(fit)
}
