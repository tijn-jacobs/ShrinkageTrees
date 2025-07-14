#' @useDynLib ShrinkageTrees, .registration = TRUE
#' @importFrom stats sd qchisq qnorm runif coef
#' @export
CausalHorseForest <- function(y,
                             X_control_train,
                             X_treat_train,
                             treatment_indicator_train,
                             X_control_test = NULL,
                             X_treat_test = NULL,
                             treatment_indicator_test = NULL,
                             status = NULL,
                             number_of_trees_control = 200,
                             number_of_trees_treat = 200,
                             N_post = 5000,
                             N_burn = 5000,
                             delayed_proposal = 5,
                             power = 2.0,
                             base = 0.95,
                             p_grow = 0.4,
                             p_prune = 0.4,
                             nu = 3,
                             q = 0.90,
                             alpha_local_control,
                             alpha_local_treat,
                             alpha_global_control = NULL,
                             alpha_global_treat = NULL,
                             tau_control = NULL,
                             tau_treat = NULL,
                             forest_wide_shrinkage_control = FALSE,
                             forest_wide_shrinkage_treat = FALSE,
                             sigma = NULL,
                             omega_control = NULL,
                             omega_treat = NULL,
                             y_mean = NULL,
                             store_parameters = FALSE,
                             store_posterior_sample_control = FALSE,
                             store_posterior_sample_treat = FALSE,
                             seed = NULL,
                             scale = "time",
                             verbose = TRUE) {
   

  # Check inputs
  if (is.null(tau_control) && is.null(alpha_global_control)) {
    stop("Specify at least one of 'tau_control' or 'alpha_global_control'.")
  }
  if (is.null(tau_treat) && is.null(alpha_global_treat)) {
    stop("Specify at least one of 'tau_treat' or 'alpha_global_treat'.")
  }

  n_train <- nrow(X_control_train)
  p_control <- ncol(X_control_train)
  p_treat <- ncol(X_treat_train)
  
  # Test data
  if (!is.null(X_control_test) && !is.null(X_treat_test)) {
    n_test <- nrow(X_control_test)
    if (ncol(X_control_test) != p_control || ncol(X_treat_test) != p_treat) {
      stop("Number of columns in X_control_test or X_treat_test does not match training data.")
    }
    X_control_test <- as.numeric(t(X_control_test))
    X_treat_test <- as.numeric(t(X_treat_test))
  } else {
    n_test <- 0
    X_control_test <- numeric(0)
    X_treat_test <- numeric(0)
    treatment_indicator_test <- integer(0)
  }
  
  # Sigma handling
  if (is.null(sigma)) {
    sigma <- sd(y)
    sigma_known <- FALSE
  } else {
    sigma_known <- TRUE
  }
  
  qchi <- qchisq(1.0 - q, nu)
  lambda <- (sigma * sigma * qchi) / nu
  
  # Seed
  if (is.null(seed)) seed <- round(runif(1, 0, 100))
  
  # Treatment indicators as integer
  treatment_indicator_train <- as.integer(treatment_indicator_train)
  if (n_test > 0) {
    treatment_indicator_test <- as.integer(treatment_indicator_test)
  }
  
  # Center or log-transform outcome
  survival <- !is.null(status)
  if (is.null(y_mean)) {
    if (survival) {
      if (scale == "time") y <- log(y)
      y_mean <- mean(y)
      y <- y - y_mean
    } else {
      y_mean <- mean(y)
      y <- y - y_mean
    }
  }

  # Scale forest with observed variance, note that this is underestimated for survival data
  if(is.null(omega_control)) {
    omega_control <- sd(y)
  }
  if(is.null(omega_treat)) {
    omega_treat <- sd(y)
  }

  # Determine prior type and param2 for control
  if (!is.null(tau_control)) {
    prior_type_control <- "fixed"
    param2_control <- tau_control
  } else if (forest_wide_shrinkage_control) {
    prior_type_control <- "horseshoe_fw"  
    param2_control <- alpha_global_control
  } else {
    prior_type_control <- "horseshoe"  
    param2_control <- alpha_global_control
  }

  # Determine prior type and param2 for treat
  if (!is.null(tau_treat)) {
    prior_type_treat <- "fixed"
    param2_treat <- tau_treat
  } else if (forest_wide_shrinkage_treat) {
    prior_type_treat <- "horseshoe_fw"  
    param2_treat <- alpha_global_treat
  } else {
    prior_type_treat <- "horseshoe"  
    param2_treat <- alpha_global_treat
  }

  
  # Call the C++ function
  causal_fit <- CausalHorseForest_cpp(
    nSEXP = n_train,
    p_treatSEXP = p_treat,
    p_controlSEXP = p_control,
    X_treatSEXP = as.numeric(t(X_treat_train)),
    X_controlSEXP = as.numeric(t(X_control_train)),
    ySEXP = y,
    status_indicatorSEXP = if (survival) status else rep(1, n_train),
    is_survivalSEXP = survival,
    treatment_indicatorSEXP = treatment_indicator_train,
    n_testSEXP = n_test,
    X_control_testSEXP = X_control_test,
    X_treat_testSEXP = X_treat_test,
    treatment_indicator_testSEXP = treatment_indicator_test,
    no_trees_treatSEXP = number_of_trees_treat,
    power_treatSEXP = power,
    base_treatSEXP = base,
    p_grow_treatSEXP = p_grow,
    p_prune_treatSEXP = p_prune,
    omega_treatSEXP = omega_treat,
    prior_type_treatSEXP = prior_type_treat,
    param1_treatSEXP = alpha_local_treat,
    param2_treatSEXP = param2_treat,
    reversible_treatSEXP = TRUE,
    eta_treatSEXP = 1,
    no_trees_controlSEXP = number_of_trees_control,
    power_controlSEXP = power,
    base_controlSEXP = base,
    p_grow_controlSEXP = p_grow,
    p_prune_controlSEXP = p_prune,
    omega_controlSEXP = omega_control,
    prior_type_controlSEXP = prior_type_control,
    param1_controlSEXP = alpha_local_control,
    param2_controlSEXP = param2_control,
    reversible_controlSEXP = TRUE,
    eta_controlSEXP = 1,
    sigma_knownSEXP = sigma_known,
    sigmaSEXP = sigma,
    lambdaSEXP = lambda,
    nuSEXP = nu,
    N_postSEXP = N_post,
    N_burnSEXP = N_burn,
    delayed_proposalSEXP = delayed_proposal,
    store_parametersSEXP = store_parameters,
    max_stored_leafsSEXP = 25,
    store_posterior_sample_controlSEXP = store_posterior_sample_control,
    store_posterior_sample_treatSEXP = store_posterior_sample_treat,
    n1SEXP = seed,
    n2SEXP = 420,
    verboseSEXP = verbose
  )
  
  # Re-transform if needed
  # Re-transform if needed
  if (survival) {
    
    if (scale == "time") {
      
      # Total
      causal_fit$train_predictions <- exp(causal_fit$train_predictions + y_mean)
      causal_fit$test_predictions <- exp(causal_fit$test_predictions + y_mean)
      
      # Control
      causal_fit$train_predictions_control <- exp(causal_fit$train_predictions_control + y_mean)
      causal_fit$test_predictions_control <- exp(causal_fit$test_predictions_control + y_mean)
      
      # Treatment
      causal_fit$train_predictions_treat <- exp(causal_fit$train_predictions_treat)
      causal_fit$test_predictions_treat <- exp(causal_fit$test_predictions_treat)
      
      # Posterior samples
      if (store_posterior_sample_control) { #
        causal_fit$train_predictions_sample_control <- exp(causal_fit$train_predictions_sample_control + y_mean)
        causal_fit$test_predictions_sample_control <- exp(causal_fit$test_predictions_sample_control + y_mean)
      }

      if (store_posterior_sample_treat) { 
        causal_fit$train_predictions_sample_treat <- exp(causal_fit$train_predictions_sample_treat)
        causal_fit$test_predictions_sample_treat <- exp(causal_fit$test_predictions_sample_treat)
      }
    } else {
      # Total
      causal_fit$train_predictions <- causal_fit$train_predictions + y_mean
      causal_fit$test_predictions <- causal_fit$test_predictions + y_mean
      
      # Control
      causal_fit$train_predictions_control <- causal_fit$train_predictions_control + y_mean
      causal_fit$test_predictions_control <- causal_fit$test_predictions_control + y_mean
      
      # Treatment
      causal_fit$train_predictions_treat <- causal_fit$train_predictions_treat
      causal_fit$test_predictions_treat <- causal_fit$test_predictions_treat
      
      # Posterior samples
      if (store_posterior_sample_control) {
        causal_fit$train_predictions_sample_control <- causal_fit$train_predictions_sample_control + y_mean
        causal_fit$test_predictions_sample_control <- causal_fit$test_predictions_sample_control + y_mean
      }

      if (store_posterior_sample_treat) {
        causal_fit$train_predictions_sample_treat <- causal_fit$train_predictions_sample_treat
        causal_fit$test_predictions_sample_treat <- causal_fit$test_predictions_sample_treat
      }
    }

    
  } else {
    # Total
    causal_fit$train_predictions <- causal_fit$train_predictions + y_mean
    causal_fit$test_predictions <- causal_fit$test_predictions + y_mean
    
    # Control
    causal_fit$train_predictions_control <- causal_fit$train_predictions_control + y_mean
    causal_fit$test_predictions_control <- causal_fit$test_predictions_control + y_mean
    
    # Treatment
    causal_fit$train_predictions_treat <- causal_fit$train_predictions_treat
    causal_fit$test_predictions_treat <- causal_fit$test_predictions_treat
    
    # Posterior samples
    if (!is.null(causal_fit$train_predictions_sample)) {
      causal_fit$train_predictions_sample <- causal_fit$train_predictions_sample + y_mean
      causal_fit$test_predictions_sample <- causal_fit$test_predictions_sample + y_mean
      causal_fit$train_predictions_sample_control <- causal_fit$train_predictions_sample_control + y_mean
      causal_fit$test_predictions_sample_control <- causal_fit$test_predictions_sample_control + y_mean
      causal_fit$train_predictions_sample_treat <- causal_fit$train_predictions_sample_treat
      causal_fit$test_predictions_sample_treat <- causal_fit$test_predictions_sample_treat
    }
  }
  
  
  return(causal_fit)
}
