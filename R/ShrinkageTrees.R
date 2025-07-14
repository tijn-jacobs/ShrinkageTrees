#' @useDynLib ShrinkageTrees, .registration = TRUE
#' @importFrom stats sd qchisq qnorm runif coef
#' @export
ShrinkageTrees <- function(y,
                      X_train,
                      X_test,
                      status = NULL,
                      number_of_trees = 200,
                      N_post = 5000,
                      N_burn = 5000,
                      delayed_proposal = 5,
                      power = 2.0,
                      base = 0.95,
                      p_grow = 0.4,
                      p_prune = 0.4,
                      nu = 3, # hyperparamater for error distribution
                      q = 0.90, # hyperparamater for error distribution
                      alpha_local, # hyperparameters for horseshoe
                      alpha_global, # hyperparameters for horseshoe
                      tau = NULL,
                      forest_wide_shrinkage = FALSE,
                      sigma = NULL,
                      binary_outcome = FALSE,
                      omega = 1,
                      store_parameters = FALSE,
                      store_posterior_sample = FALSE,
                      seed = NULL,
                      scale = "time",
                      verbose = TRUE
                      ) { 
  
  # Retrieve size information of the training- and test data.
  n_train <- nrow(X_train)         # number of observations in training data
  p_features <- ncol(X_train)      # number of features
  n_test <- nrow(X_train)           # number of test observations
  

  # Check if X_test is provided, and if it has the correct dimensions
  if (!is.null(X_test)) {
    n_test <- nrow(X_test)         # number of test observations
    
    # Ensure X_test has the correct number of covariates
    if (ncol(X_test) != p_features) {
      stop("The number of covariates in X_test must match the number in
            X_train.")
    }
    
    # Flatten the test data (should be contiguous)
    X_test <- as.numeric(t(X_test))
  } else {
    
    # Fill up empty test data by the mean of the training data 
    n_test <- 1                    
    X_test <- as.numeric(colMeans(X_train))
  }
  # Compute empirical prior parameters for error distribution  
  if(is.null(sigma)) {
    sigma <- sd(y)
    sigma_known = FALSE
  } else {
    sigma_known = TRUE
  }
  qchi = qchisq(1.0 - q, nu)
  lambda = (sigma*sigma*qchi)/nu # lambda parameter for sigma prior
  
  # Ensure data is in the correct format for C++
  N_post <- as.integer(N_post)[1]
  N_burn <- as.integer(N_burn)[1]
  power <- as.numeric(power)[1]
  base <- as.numeric(base)[1]
  p_grow <- as.numeric(p_grow)[1]
  p_prune <- as.numeric(p_prune)[1]
  
  
  # Flatten the training and test data (should be continguous)
  X_train <- as.numeric(t(X_train))
  X_test <- as.numeric(t(X_test))
  
  # Set a random seed
  if(is.null(seed)) seed <- round(runif(1, 0, 100))

  if (!is.null(status)) { # Outcome is survival

    # Log transform the data
    y <- as.numeric(y)
    if (scale == "time") y <- log(y)
    y_mean = mean(y)
    y <- y - y_mean
    survival <- TRUE

    if(is.null(tau)) {
      
      # No specified tau means a random global shrinkage parameter.
      # This can be on the tree level or forest wide.
      if(forest_wide_shrinkage) {
        fit <- HorseTrees_cpp(nSEXP = n_train,
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
                          sigmaSEXP = sigma,
                          sigma_knownSEXP = sigma_known,
                          omegaSEXP = omega,
                          param1SEXP = alpha_local,
                          param2SEXP = alpha_global,
                          prior_typeSEXP = "horseshoe_fw", # Indicates Horseshoe_fw ScaleMixture
                          reversibleSEXP = TRUE,
                          store_parametersSEXP = store_parameters,
                          store_posterior_sampleSEXP = store_posterior_sample,
                          n1SEXP = seed, 
                          n2SEXP = 420,
                          verboseSEXP = verbose)
      } else {
        fit <- HorseTrees_cpp(nSEXP = n_train,
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
                          sigmaSEXP = sigma,
                          sigma_knownSEXP = sigma_known,
                          omegaSEXP = omega,
                          param1SEXP = alpha_local,
                          param2SEXP = alpha_global,
                          prior_typeSEXP = "horseshoe", # Indicates Horseshoe ScaleMixture
                          reversibleSEXP = TRUE,
                          store_parametersSEXP = store_parameters,
                          store_posterior_sampleSEXP = store_posterior_sample,
                          n1SEXP = seed, 
                          n2SEXP = 420,
                          verboseSEXP = verbose)
      }
    #} else { 
      if (scale == "time") {

        # Re-transform the data
        fit$train_predictions = exp(fit$train_predictions + y_mean)
        fit$test_predictions = exp(fit$test_predictions + y_mean)
        fit$train_predictions_sample = exp(fit$train_predictions_sample + y_mean)
        fit$test_predictions_sample = exp(fit$test_predictions_sample + y_mean)
      } else {
        # Re-transform the data
        fit$train_predictions = fit$train_predictions + y_mean
        fit$test_predictions = fit$test_predictions + y_mean
        fit$train_predictions_sample = fit$train_predictions_sample + y_mean
        fit$test_predictions_sample = fit$test_predictions_sample + y_mean
      }


    } else {
      
      # tau is specified; fit a Horseshoe with fixed global shrinkage
      fit <- HorseTrees_cpp(nSEXP = n_train,
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
                        sigmaSEXP = sigma,
                        sigma_knownSEXP = sigma_known,
                        omegaSEXP = omega,
                        param1SEXP = alpha_local,
                        param2SEXP = tau,
                        prior_typeSEXP = "halfcauchy", # Indicates (or Horseshoe with fixed global shrink. param.)
                        reversibleSEXP = TRUE,
                        store_parametersSEXP = store_parameters,
                        store_posterior_sampleSEXP = store_posterior_sample,
                        n1SEXP = seed,        
                        n2SEXP = 420,
                        verboseSEXP = verbose          
      )
    }

  } else { # Outcome is binary or continuous

    status <- rep(1, n_train) # unused
    survival <- FALSE
    
    if (binary_outcome) {
          # Only center(!) the outcomes for continuous data
    y <- as.numeric(y)
    latent_threshold = qnorm(mean(y))

    fit <- probitHorseTrees_cpp(nSEXP = n_train,
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
                        omegaSEXP = omega,
                        latent_thresholdSEXP = latent_threshold,
                        param1SEXP = alpha_local,
                        param2SEXP = alpha_global,
                        prior_typeSEXP = "horseshoe",
                        reversibleSEXP = TRUE,
                        store_posterior_sampleSEXP = store_posterior_sample,
                        n1SEXP = seed, 
                        n2SEXP = 420,
                        verboseSEXP = verbose 
      )

  } else {

      # Only center(!) the outcomes for continuous data
      y <- as.numeric(y)
      y_mean = mean(y)
      y_sd = sd(y)
      y <- y - y_mean

    if(is.null(tau)) {
      
      # No specified tau means a random global shrinkage parameter.
      # This can be on the tree level or forest wide.
      if(forest_wide_shrinkage) {
        fit <- HorseTrees_cpp(nSEXP = n_train,
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
                          sigmaSEXP = sigma,
                          sigma_knownSEXP = sigma_known,
                          omegaSEXP = omega,
                          param1SEXP = alpha_local,
                          param2SEXP = alpha_global,
                          prior_typeSEXP = "horseshoe_fw", # Indicates Horseshoe_fw ScaleMixture
                          reversibleSEXP = TRUE,
                          store_parametersSEXP = store_parameters,
                          store_posterior_sampleSEXP = store_posterior_sample,
                          n1SEXP = seed, 
                          n2SEXP = 420,
                          verboseSEXP = verbose)
      } else {
        fit <- HorseTrees_cpp(nSEXP = n_train,
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
                          sigmaSEXP = sigma,
                          sigma_knownSEXP = sigma_known,
                          omegaSEXP = omega,
                          param1SEXP = alpha_local,
                          param2SEXP = alpha_global,
                          prior_typeSEXP = "horseshoe", # Indicates Horseshoe ScaleMixture
                          reversibleSEXP = TRUE,
                          store_parametersSEXP = store_parameters,
                          store_posterior_sampleSEXP = store_posterior_sample,
                          n1SEXP = seed, 
                          n2SEXP = 420,
                          verboseSEXP = verbose)
      }


    } else {
      
      # tau is specified; fit a Horseshoe with fixed global shrinkage
      fit <- HorseTrees_cpp(nSEXP = n_train,
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
                        sigmaSEXP = sigma,
                        sigma_knownSEXP = sigma_known,
                        omegaSEXP = omega,
                        param1SEXP = alpha_local,
                        param2SEXP = tau,
                        prior_typeSEXP = "halfcauchy", # Indicates (or Horseshoe with fixed global shrink. param.)
                        reversibleSEXP = TRUE,
                        store_parametersSEXP = store_parameters,
                        store_posterior_sampleSEXP = store_posterior_sample,
                        n1SEXP = seed,        
                        n2SEXP = 420,
                        verboseSEXP = verbose)
    }

    # Recenter the data
    fit$train_predictions = fit$train_predictions + y_mean
    fit$test_predictions = fit$test_predictions + y_mean
    fit$train_predictions_sample = fit$train_predictions_sample + y_mean
    fit$test_predictions_sample = fit$test_predictions_sample + y_mean

  }
  }
  

  
  return(fit)
}