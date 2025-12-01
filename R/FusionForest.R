#' Fusion Forests
#'
#' TBD
#'
#' @importFrom Rcpp evalCpp
#' @useDynLib ShrinkageTrees, .registration = TRUE
#' @importFrom stats sd qchisq qnorm runif coef
#' @export
FusionShrinkageForest <- function(y,
                                  status = NULL,
                                  X_train_control,
                                  X_train_treat,
                                  treatment_indicator_train,
                                  source_indicator_train,
                                  X_test_control = NULL,
                                  X_test_treat = NULL,
                                  X_test_deconf = NULL,
                                  treatment_indicator_test = NULL,
                                  source_indicator_test = NULL,
                                  outcome_type = "continuous",
                                  timescale = "time",
                                  number_of_trees_control = 200,
                                  number_of_trees_treat = 200,
                                  number_of_trees_deconf = 200,
                                  prior_type_control = "horseshoe",
                                  prior_type_treat   = "horseshoe",
                                  prior_type_deconf  = "horseshoe",
                                  local_hp_control,
                                  local_hp_treat,
                                  local_hp_deconf,
                                  global_hp_control = NULL,
                                  global_hp_treat   = NULL,
                                  global_hp_deconf  = NULL,
                                  power = 2.0,
                                  base = 0.95,
                                  p_grow = 0.4,
                                  p_prune = 0.4,
                                  nu = 3,
                                  eta_commensurate = 0.0,
                                  q = 0.90,
                                  sigma = NULL,
                                  N_post = 5000,
                                  N_burn = 5000,
                                  delayed_proposal = 5,
                                  store_posterior_sample = FALSE,
                                  seed = NULL,
                                  verbose = TRUE
) {


  
  # Check prior_type value
  allowed_prior <- c("horseshoe", "horseshoe_fw", "horseshoe_EB", "half-cauchy", "standard")
  if (!prior_type_control %in% allowed_prior) {
    stop("Invalid prior_type_control Choose 'horseshoe', 'horseshoe_fw', 
         'horseshoe_EB', 'half-cauchy' or 'standard'.")
  }
  
  if (!prior_type_treat %in% allowed_prior) {
    stop("Invalid prior_type_treat. Choose 'horseshoe', 'horseshoe_fw', 
         'horseshoe_EB', 'half-cauchy'or 'standard'.")
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
  
  if (prior_type_control %in% c("half-cauchy", "standard")) {
    if (is.null(local_hp_control)) {
      stop("For prior_type = 'half-cauchy' or 'standard', you must provide local_hp_control.")
    }
    if (!is.null(global_hp_control)) {
      warning("global_hp_control is ignored for 'half-cauchy' or 'standard'. If you want to 
              fix the global parameter, consider using 'horseshoe_EB'.")
    }
    global_hp_control <- 1
    prior_type_control <- "halfcauchy"
  }
  
  if (prior_type_treat %in% c("half-cauchy", "standard")) {
    if (is.null(local_hp_treat)) {
      stop("For prior_type_treat = 'half-cauchy' or 'standard', you must provide 
           local_hp_treat.")
    }
    if (!is.null(global_hp_treat)) {
      warning("global_hp_treat is ignored for 'half-cauchy' or 'standard'. If you want to fix 
              the global parameter, consider using 'horseshoe_EB'.")
    }
    global_hp_treat <- 1
    prior_type_treat <- "halfcauchy"
  }

  if (prior_type_deconf %in% c("half-cauchy", "standard")) {
    if (is.null(local_hp_deconf)) {
      stop("For prior_type_deconf = 'half-cauchy' or 'standard', you must provide 
           local_hp_deconf.")
    }
    if (!is.null(global_hp_deconf)) {
      warning("global_hp_deconf is ignored for 'half-cauchy' or 'standard'. If you want to fix 
              the global parameter, consider using 'horseshoe_EB'.")
    }
    global_hp_deconf <- 1
    prior_type_deconf <- "halfcauchy"
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
  
  
  # # Retrieve dimensions of training data
  # n_train <- nrow(X_train_control)

  # p_control <- ncol(X_train_control)
  # p_treat   <- ncol(X_train_treat)
  
  # # Check matching row numbers
  # if (nrow(X_train_control) != length(y)) {
  #   stop("X_train_control rows must match length of y.")
  # }
  # if (nrow(X_train_treat) != length(y)) {
  #   stop("X_train_treat rows must match length of y.")
  # }
  # if (length(treatment_indicator_train) != length(y)) {
  #   stop("treatment_indicator_train must match length of y.")
  # }
  # if (length(source_indicator_train) != length(y)) {
  #   stop("source_indicator_train must match length of y.")
  # }

  # # If test provided, check
  # if (!is.null(X_test_control) && !is.null(X_test_treat)) {
  #   n_test <- nrow(X_test_control)
  #   if (ncol(X_test_control) != p_control || ncol(X_test_treat) != p_treat) {
  #     stop("Number of columns in X_test_control or X_test_treat does not match 
  #          training data.")
  #   }
  #   if (!is.null(treatment_indicator_test) && 
  #       length(treatment_indicator_test) != n_test) {
  #     stop("treatment_indicator_test length must match number of test rows.")
  #   }
  #   X_test_control <- as.numeric(t(X_test_control))
  #   X_test_treat <- as.numeric(t(X_test_treat))
  # } else {
  #   n_test <- 1
  #   X_test_control <- as.numeric(colMeans(X_train_control))
  #   X_test_treat <- as.numeric(colMeans(X_train_treat))
  #   X_test_deconf  <- X_test_control                    # <- default for deconf
  #   treatment_indicator_test <- as.integer(1L)
  #   source_indicator_test    <- as.integer(1L)          # all RCT by default
  # }

  # ## CHECK LATER!

  # # source indicator for test: default to all RCT if missing
  # if (is.null(source_indicator_test)) {
  #   source_indicator_test <- rep.int(1L, n_test)
  # } else {
  #   source_indicator_test <- as.integer(source_indicator_test)
  #   if (length(source_indicator_test) != n_test)
  #     stop("source_indicator_test length must match number of test rows.")
  #   if (!all(source_indicator_test %in% c(0L, 1L)))
  #     stop("source_indicator_test must be 0 (OS) or 1 (RCT).")
  # }

  # # Treatment indicators as integer
  # treatment_indicator_train <- as.integer(treatment_indicator_train)
  # if (is.null(treatment_indicator_test)) {
  #   treatment_indicator_test <- as.integer(rep(1, n_test))
  # } else {
  #   treatment_indicator_test <- as.integer(treatment_indicator_test)
  # }

  # # source indicators as integer
  # if (!all(source_indicator_train %in% c(0L,1L))) stop("source_indicator_train must be 0 (OS) or 1 (RCT).")
  # source_indicator_train <- as.integer(source_indicator_train)
  # n_deconf <- sum(source_indicator_train == 0L)
  # if (n_deconf <= 0) stop("n_deconf == 0 not allowed (need OS rows).")
  # X_train_deconf <- X_train_control[source_indicator_train == 0L, , drop = FALSE]
  # p_deconf <- ncol(X_train_deconf)
  # X_train_deconf <- as.numeric(t(X_train_deconf))

  #   # deconf test matrix: default to X_test_control if missing
  # if (is.null(X_test_deconf)) {
  #   X_test_deconf <- X_test_control
  # } else {
  #   if (ncol(X_test_deconf) != p_deconf)
  #     stop("X_test_deconf must have ", p_deconf, " columns.")
  #   X_test_deconf <- as.numeric(t(X_test_deconf))
  # }

  
  # # Force data types to be numeric and plain arrays
  # N_post <- as.integer(N_post)[1]
  # N_burn <- as.integer(N_burn)[1]
  # power <- as.numeric(power)[1]
  # base <- as.numeric(base)[1]
  # p_grow <- as.numeric(p_grow)[1]
  # p_prune <- as.numeric(p_prune)[1]
  # X_train_treat <- as.numeric(t(X_train_treat))
  # X_train_control <- as.numeric(t(X_train_control))

  ## TEST FROM HERE ##

  ## ---------------------------
  ## Data preparation & checks
  ## ---------------------------

  # Basic dimensions
  n_train   <- nrow(X_train_control)
  p_control <- ncol(X_train_control)
  p_treat   <- ncol(X_train_treat)

  # Consistency of training sizes
  if (nrow(X_train_control) != length(y))
    stop("X_train_control rows must match length of y.")
  if (nrow(X_train_treat)   != length(y))
    stop("X_train_treat rows must match length of y.")
  if (length(treatment_indicator_train) != length(y))
    stop("treatment_indicator_train must match length of y.")
  if (length(source_indicator_train)    != length(y))
    stop("source_indicator_train must match length of y.")

  # Coerce training indicators to integer and validate source {0,1}
  treatment_indicator_train <- as.integer(treatment_indicator_train)
  if (!all(source_indicator_train %in% c(0L, 1L)))
    stop("source_indicator_train must be 0 (OS) or 1 (RCT).")
  source_indicator_train <- as.integer(source_indicator_train)

  # OS subset for deconf forest
  n_deconf <- sum(source_indicator_train == 0L)
  if (n_deconf <= 0L)
    stop("n_deconf == 0 not allowed (need at least one OS row).")
  X_train_deconf <- X_train_control[source_indicator_train == 0L, , drop = FALSE]
  p_deconf       <- ncol(X_train_deconf)

  # --- Test data handling (validate as matrices first, then flatten) ---

  if (!is.null(X_test_control) && !is.null(X_test_treat)) {
    # Ensure matrix inputs
    if (!is.matrix(X_test_control) || !is.matrix(X_test_treat))
      stop("X_test_control and X_test_treat must be matrices when provided.")
    n_test <- nrow(X_test_control)
    if (nrow(X_test_treat) != n_test)
      stop("X_test_control and X_test_treat must have the same number of rows.")
    if (ncol(X_test_control) != p_control || ncol(X_test_treat) != p_treat)
      stop("X_test_* column counts must match training: p_control and p_treat.")

    # Default/validate test source indicators
    if (is.null(source_indicator_test)) {
      source_indicator_test <- rep.int(1L, n_test)  # default: all RCT
    } else {
      source_indicator_test <- as.integer(source_indicator_test)
      if (length(source_indicator_test) != n_test)
        stop("source_indicator_test length must match number of test rows.")
      if (!all(source_indicator_test %in% c(0L, 1L)))
        stop("source_indicator_test must be 0 (OS) or 1 (RCT).")
    }

    # Default/validate test treatment indicators
    if (is.null(treatment_indicator_test)) {
      treatment_indicator_test <- rep.int(1L, n_test)
    } else {
      treatment_indicator_test <- as.integer(treatment_indicator_test)
      if (length(treatment_indicator_test) != n_test)
        stop("treatment_indicator_test length must match number of test rows.")
      if (!all(treatment_indicator_test %in% c(0L, 1L)))
        stop("treatment_indicator_test must be 0/1.")
    }

    # Deconf test block: default to control; else validate shape
    if (is.null(X_test_deconf)) {
      X_test_deconf <- X_test_control
    } else {
      if (!is.matrix(X_test_deconf))
        stop("X_test_deconf must be a matrix when provided.")
      if (nrow(X_test_deconf) != n_test)
        stop("X_test_deconf must have the same number of rows as X_test_control.")
      if (ncol(X_test_deconf) != p_deconf)
        stop("X_test_deconf must have ", p_deconf, " columns.")
    }

    # Flatten test matrices after checks (C++ expects column-major vector via t(.))
    X_test_control <- as.numeric(t(X_test_control))
    X_test_treat   <- as.numeric(t(X_test_treat))
    X_test_deconf  <- as.numeric(t(X_test_deconf))

  } else {
    # No explicit test data: use single-row defaults at train means
    n_test            <- 1L
    X_test_control    <- matrix(colMeans(X_train_control), nrow = 1L)
    X_test_treat      <- matrix(colMeans(X_train_treat),   nrow = 1L)
    X_test_deconf     <- X_test_control
    treatment_indicator_test <- 1L
    source_indicator_test    <- 1L

    # Flatten the single-row test matrices
    X_test_control <- as.numeric(t(X_test_control))
    X_test_treat   <- as.numeric(t(X_test_treat))
    X_test_deconf  <- as.numeric(t(X_test_deconf))
  }

  # Flatten training matrices after all row/column checks
  X_train_control  <- as.numeric(t(X_train_control))
  X_train_treat    <- as.numeric(t(X_train_treat))
  X_train_deconf   <- as.numeric(t(X_train_deconf))

  # Scalarize tunables to base types C++ expects
  N_post  <- as.integer(N_post)[1L]
  N_burn  <- as.integer(N_burn)[1L]
  power   <- as.numeric(power)[1L]
  base    <- as.numeric(base)[1L]
  p_grow  <- as.numeric(p_grow)[1L]
  p_prune <- as.numeric(p_prune)[1L]


  ## TEST TO HERE ##

  
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
    fit <- FusionForest_cpp(
      nSEXP = n_train,
      p_treatSEXP = p_treat,
      p_controlSEXP = p_control,
      X_train_treatSEXP = X_train_treat,
      X_train_controlSEXP = X_train_control,
      ySEXP = y,
      status_indicatorSEXP = status,
      is_survivalSEXP = survival,
      treatment_indicatorSEXP = treatment_indicator_train,
      source_indicatorSEXP = source_indicator_train,
      n_testSEXP = n_test,
      X_test_controlSEXP = X_test_control,
      X_test_treatSEXP = X_test_treat,
      X_test_deconfSEXP = X_test_deconf,
      treatment_indicator_testSEXP = treatment_indicator_test,
      source_indicator_testSEXP = source_indicator_test,
      n_deconfSEXP = n_deconf,
      p_deconfSEXP = p_deconf,
      X_train_deconfSEXP = X_train_deconf,
      no_trees_deconfSEXP = number_of_trees_deconf,
      power_deconfSEXP = power,
      base_deconfSEXP  = base,
      p_grow_deconfSEXP = p_grow,
      p_prune_deconfSEXP = p_prune,
      omega_deconfSEXP = 1/2,
      prior_type_deconfSEXP = prior_type_deconf,
      param1_deconfSEXP = local_hp_deconf,
      param2_deconfSEXP = global_hp_deconf,
      reversible_deconfSEXP = TRUE,
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
      eta_commensurateSEXP = eta_commensurate,
      N_postSEXP = N_post,
      N_burnSEXP = N_burn,
      delayed_proposalSEXP = delayed_proposal,
      store_parametersSEXP = FALSE,
      max_stored_leavesSEXP = 1,
      store_posterior_sampleSEXP = store_posterior_sample,
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
      
      # Deconfounding
      fit$train_predictions_deconf <- exp(fit$train_predictions_deconf * sigma_hat)
      fit$test_predictions_deconf  <- exp(fit$test_predictions_deconf  * sigma_hat)

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
        fit$train_predictions_sample_deconf <- 
          exp(fit$train_predictions_sample_deconf * sigma_hat)
        fit$test_predictions_sample_deconf  <- 
          exp(fit$test_predictions_sample_deconf  * sigma_hat)
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

      # Deconfounding
      fit$train_predictions_deconf <- fit$train_predictions_deconf * sigma_hat
      fit$test_predictions_deconf  <- fit$test_predictions_deconf  * sigma_hat
      
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
    fit <- FusionForest_cpp(
      nSEXP = n_train,
      p_treatSEXP = p_treat,
      p_controlSEXP = p_control,
      X_train_treatSEXP = X_train_treat,
      X_train_controlSEXP = X_train_control,
      ySEXP = y,
      status_indicatorSEXP = status,         # dummy 1s here (continuous)
      is_survivalSEXP = FALSE,
      treatment_indicatorSEXP = treatment_indicator_train,
      source_indicatorSEXP = source_indicator_train,   

      n_testSEXP = n_test,
      X_test_controlSEXP = X_test_control,
      X_test_treatSEXP = X_test_treat,
      X_test_deconfSEXP = X_test_deconf,               
      treatment_indicator_testSEXP = treatment_indicator_test,
      source_indicator_testSEXP = source_indicator_test, 
      n_deconfSEXP = n_deconf,                         
      p_deconfSEXP = p_deconf,                         
      X_train_deconfSEXP = X_train_deconf,      

      # --- deconf forest hyperparams ---
      no_trees_deconfSEXP = number_of_trees_deconf,
      power_deconfSEXP = power,
      base_deconfSEXP  = base,
      p_grow_deconfSEXP = p_grow,
      p_prune_deconfSEXP = p_prune,
      omega_deconfSEXP = 1/2,
      prior_type_deconfSEXP = prior_type_deconf,
      param1_deconfSEXP = local_hp_deconf,
      param2_deconfSEXP = global_hp_deconf,
      reversible_deconfSEXP = TRUE,

      # --- treatment forest hyperparams ---
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

      # --- control (prognostic) forest hyperparams ---
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

      # --- noise + borrowing params ---
      sigma_knownSEXP = sigma_known,
      sigmaSEXP = sigma_hat,
      lambdaSEXP = lambda,
      nuSEXP = nu,     
      eta_commensurateSEXP = eta_commensurate,    

      # --- MCMC / storage / misc ---
      N_postSEXP = N_post,
      N_burnSEXP = N_burn,
      delayed_proposalSEXP = delayed_proposal,
      store_parametersSEXP = FALSE,
      max_stored_leavesSEXP = 1,
      store_posterior_sampleSEXP = store_posterior_sample,
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

    # Deconfounding
    fit$train_predictions_deconf <- fit$train_predictions_deconf * sigma_hat
    fit$test_predictions_deconf  <- fit$test_predictions_deconf  * sigma_hat
    
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
      fit$train_predictions_sample_deconf <- 
        fit$train_predictions_sample_deconf * sigma_hat
      fit$test_predictions_sample_deconf  <- 
        fit$test_predictions_sample_deconf  * sigma_hat
    }
    
  }
  
  # remove burn-in of sigma
  if (!sigma_known) {
    fit$sigma <- fit$sigma[-(1:N_burn)]
  } 

  return(fit)
}