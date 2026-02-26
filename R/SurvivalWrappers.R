#' SurvivalBART 
#'
#' Fits an Accelerated Failure Time (AFT) model using the original BART prior:
#' We fit a BART model to f:
#'  log(y) = f(x) + ε.
#'
#' @param y Outcome vector of right-censored (non-negative) survival times.
#' @param status Event indicator (1 = event, 0 = censored). Required for survival.
#' @param X_train Design matrix for training.
#' @param X_test Optional test matrix (defaults to column means of X_train).
#' @param number_of_trees Number of trees (default 200).
#' @param k Hyperparameter for the prior variance of the step heights.
#' @param N_post Number of posterior samples to draw.
#' @param N_burn Number of burn-in iterations.
#' @param verbose Print sampling progress.
#'
#' @details
#' This function provides a simplified interface for fitting an AFT model with 
#' the classical BART prior using the internal \code{ShrinkageTrees} engine.
#' Right-censored survival data are handled via data augmentation within the 
#' tree sampler.
#'
#' Compared to the full \code{ShrinkageTrees} function, \code{SurvivalBART} fixes several hyperparameters.
#' This wrapper is intended as a convenient survival-specific version of BART for 
#' users who do not need advanced shrinkage priors.
#' @references 
#' Chipman, H. A., George, E. I., & McCulloch, R. E. (2010).  
#' *BART: Bayesian Additive Regression Trees.* Annals of Applied Statistics.
#'  
#' @export
SurvivalBART <- function(
  time, 
  status, 
  X_train, 
  X_test = NULL,
  timescale = "time",
  number_of_trees = 200,
  k = 2.0,
  N_post = 1000, 
  N_burn = 1000,
  verbose = TRUE
) {

  if (timescale == "time") {
    local_hp <- (max(log(time)) - min(log(time))) / (k * 2 * sqrt(number_of_trees))
  } else {
    local_hp <- (max(time) - min(time)) / (k * 2 * sqrt(number_of_trees))
  }
  ShrinkageTrees(
    y = time,
    status = status,
    X_train = X_train,
    X_test = X_test,
    outcome_type = "right-censored",
    timescale = timescale,
    number_of_trees = number_of_trees,
    prior_type = "standard",
    local_hp = local_hp,
    N_post = N_post,
    N_burn = N_burn,
    store_posterior_sample = TRUE,
    verbose = verbose
  )
}


#' SurvivalDART
#'
#' Fits an AFT model using the Dirichlet splitting prior that encourages sparse
#' feature selection via a Beta–Dirichlet hierarchy.
#'
#' @inheritParams SurvivalBART
#' @param a_dirichlet,b_dirichlet Beta hyperparameters for sparsity.
#' @param rho_dirichlet Expected number of active predictors (default = ncol(X)).
#'
#' @return A fitted AFT-DART object.
#' @export
SurvivalDART <- function(
  time, 
  status, 
  X_train, 
  X_test = NULL,
  timescale = "time",
  number_of_trees = 200,
  a_dirichlet = 0.5,
  b_dirichlet = 1.0,
  rho_dirichlet = NULL,
  k = 2.0,
  N_post = 1000,
  N_burn = 1000,
  verbose = TRUE
) {

  if (is.null(rho_dirichlet)) rho_dirichlet <- ncol(X_train)
  if (timescale == "time") {
    local_hp <- (max(log(time)) - min(log(time))) / (k * 2 * sqrt(number_of_trees))
  } else {
    local_hp <- (max(time) - min(time)) / (k * 2 * sqrt(number_of_trees))
  }
  ShrinkageTrees(
    y = time,
    status = status,
    X_train = X_train,
    X_test = X_test,
    outcome_type = "right-censored",
    timescale = timescale,
    number_of_trees = number_of_trees,
    prior_type = "dirichlet",
    local_hp = local_hp,
    a_dirichlet = a_dirichlet,
    b_dirichlet = b_dirichlet,
    rho_dirichlet = rho_dirichlet,
    N_post = N_post,
    N_burn = N_burn,
    store_posterior_sample = TRUE,
    verbose = verbose
  )
}


#' SurvivalBCF (Bayesian Causal Forest; Hahn et al. 2020)
#'
#' A simplified wrapper for fitting an AFT version of BCF:
#' Y = f(x) + a * τ(x) + ε, with separate control/treatment forests.
#'
#' @inheritParams SurvivalBART
#' @param W Treatment indicator (0/1) for training.
#' @param number_of_trees_mu Trees for control forest.
#' @param number_of_trees_tau Trees for treatment effect forest.
#' @param k_mu,k_tau Leaf prior scales.
#'
#' @export
SurvivalBCF <- function(
  time,
  status,
  X_train, 
  treatment,
  timescale = "time",
  propensity = NULL,
  N_post = 1000,
  N_burn = 1000,
  verbose = TRUE
) {
  if(!is.null(propensity)) {
    X_train_control <- cbind(propensity, X_train)
  } else {
    X_train_control <- X_train
  }
  info <- censored_info(time, status)
  sd_time <- info$sd
  CausalShrinkageForest(
    y = time,
    status = status,
    X_train_control = X_train_control,
    X_train_treat = X_train,
    treatment_indicator_train = treatment,
    outcome_type = "right-censored",
    timescale = timescale,
    prior_type_control = "standard-halfcauchy",
    prior_type_treat = "standard",
    local_hp_treat = 1 / (sd_time * sqrt(50) * 0.674),
    power_control = 2, base_control = 0.95,
    power_treat = 3, base_treat = 0.25,
    number_of_trees_control = 200,
    number_of_trees_treat = 50,
    N_post = N_post,
    N_burn = N_burn,
    store_posterior_sample = TRUE,
    verbose = verbose
  )
}

#' SurvivalShrinkageBCF (BCF + Dirichlet Sparsity)
#'
#' This model generalizes BCF by including:
#'  - Dirichlet splitting priors for sparsity (optional).
#'
#' @inheritParams SurvivalBCF
#' @param prior_mu,prior_tau Prior types for control & treatment forests.
#'
#' @export
SurvivalShrinkageBCF <- function(
  time, 
  status, 
  X_train,  
  treatment,
  timescale = "time",
  propensity = NULL, 
  a_dir = 0.5,
  b_dir = 1.0,
  N_post = 1000, 
  N_burn = 1000,
  verbose = TRUE
) {

  if (!is.null((propensity))) {
    X_control <- cbind(propensity, X_train)
    rho_dir_control <- ncol(X_control)
    rho_dir_treat <- ncol(X_train)
  } else {
    X_control <- X_train
    rho_dir_treat <- rho_dir_control <- ncol(X_train)
  }

  info <- censored_info(time, status)
  sd_time <- info$sd

  CausalShrinkageForest(
    y = time,
    status = status,
    X_train_control = X_control,
    X_train_treat = X_train,
    treatment_indicator_train = treatment,
    outcome_type = "right-censored",
    timescale = timescale,
    prior_type_control = "dirichlet-halfcauchy",
    prior_type_treat = "dirichlet",
    local_hp_treat = 1 / (sd_time * sqrt(50) * 0.674),
    power_control = 2, base_control = 0.95,
    power_treat = 3, base_treat = 0.25,
    number_of_trees_control = 200,
    number_of_trees_treat = 50,
    a_dirichlet_control = a_dir,
    b_dirichlet_control = b_dir,
    rho_dirichlet_control = rho_dir_control,
    a_dirichlet_treat = a_dir,
    b_dirichlet_treat = b_dir,
    rho_dirichlet_treat = rho_dir_treat,
    N_post = N_post,
    N_burn = N_burn,
    store_posterior_sample = TRUE,
    verbose = verbose
  )
}
