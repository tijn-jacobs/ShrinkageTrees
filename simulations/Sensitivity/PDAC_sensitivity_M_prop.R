## Simulations based on PDAC data example for reviewer 2
library(ShrinkageTrees)
library(foreach)
library(doParallel)


compute_cate_detection <- function(cate_ci, true_cate, delta) {
  idx <- abs(true_cate) > delta
  if (!any(idx)) return(NA_real_)
  detected <- cate_ci[1, idx] > 0 | cate_ci[2, idx] < 0
  mean(detected)
}

evaluate_CHF <- function(data, use_uncensored = FALSE, N_post = N_post,
                         N_burn = N_burn, sigma_known = FALSE, ...) {
  if (!requireNamespace("ShrinkageTrees", quietly = TRUE)) {
    stop("Package 'ShrinkageTrees' is required.")
  }
  
  # Determine outcome and censoring indicator
  if (use_uncensored) {
    y <- data$uncensored_event_times
    status <- NULL
  } else {
    y <- data$follow_up
    status <- data$status
  }
  
  if (sigma_known) {
    sigma <- data$obs_sigma
  } else {
    sigma <- NULL
  }
  
  X1 <- cbind(data$propensity, data$X_train)
  X2 <- matrix(colMeans(X1), nrow = 1)
  
  # Fit the CHF model
  fit <- ShrinkageTrees::CausalShrinkageForest(
    y = y,
    status = status,
    outcome_type = "right-censored",
    X_train_control = X1,
    X_train_treat = data$X_train,
    treatment_indicator_train = data$treatment,
    store_posterior_sample = TRUE,
    timescale = "log",
    verbose = FALSE,
    prior_type_control = "horseshoe",
    prior_type_treat = "horseshoe",
    sigma = sigma,
    N_post = N_post,
    N_burn = N_burn,
    ...
  )
  
  # Posterior samples of τ(x)
  cate_samples <- fit$train_predictions_sample_treat
  
  # Pointwise quantities
  cate_mean <- fit$train_predictions_treat
  cate_rpehe <- sqrt(mean((cate_mean - data$true_cate)^2))
  cate_abs_bias <- mean(abs(cate_mean - data$true_cate))
  
  # Per-observation posterior intervals for τ(x)
  cate_ci <- apply(cate_samples, 2, quantile, probs = c(0.025, 0.975))
  cate_coverage <- mean(cate_ci[1,] <= data$true_cate & cate_ci[2,] >= data$true_cate)
  cate_ci_length <- mean(cate_ci[2,] - cate_ci[1,])
  
  # ATE estimate and uncertainty
  ate_samples <- rowMeans(cate_samples)
  ate_estimate <- mean(ate_samples)
  ate_ci <- quantile(ate_samples, probs = c(0.025, 0.975))
  ate_coverage <- unname(ate_ci[1] <= data$true_ate & ate_ci[2] >= data$true_ate)
  ate_ci_length <- unname(diff(ate_ci))
  
  delta <- 0.25 * sd(data$true_cate)
  cate_detection_power <- compute_cate_detection(
    cate_ci   = cate_ci,
    true_cate = data$true_cate,
    delta     = delta
  )
  
  # RMSE for total outcome prediction
  rmse <- sqrt(mean((fit$train_predictions - data$true_event_times)^2))
  
  # Compute the C-index
  if (requireNamespace("survival", quietly = TRUE)) {
    c_index <- survival::concordance(
      survival::Surv(exp(data$follow_up), data$status) ~ exp(fit$train_predictions)
    )$concordance
  } else {
    c_index <- NA
  }
  
  return(list(
    postmean_sigma = mean(fit$sigma),
    ate = ate_estimate,
    ate_ci = ate_ci,
    ate_coverage = ate_coverage,
    ate_ci_length = ate_ci_length,
    cate_rpehe = cate_rpehe,
    cate_abs_bias  = cate_abs_bias,
    cate_coverage = cate_coverage,
    cate_ci_length = cate_ci_length,
    rmse = rmse,
    cate_detection_power  = cate_detection_power,
    c_index = c_index
  ))
}

evaluate_CHF_CV <- function(data,
                            param_grid,
                            K = 5,
                            number_of_trees = 200,
                            N_post = 1000,
                            N_burn = 1000,
                            use_uncensored = FALSE,
                            sigma_known = FALSE,
                            seed = 1) {
  set.seed(seed)
  
  if (!requireNamespace("ShrinkageTrees", quietly = TRUE)) {
    stop("Package 'ShrinkageTrees' is required.")
  }
  
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package 'survival' is required for computing C-index.")
  }
  
  # Determine outcome
  if (use_uncensored) {
    y_full <- data$uncensored_event_times
    status_full <- rep(1, length(data$status))
  } else {
    y_full <- data$follow_up
    status_full <- data$status
  }
  
  if (sigma_known) {
    sigma <- data$obs_sigma
  } else {
    sigma <- NULL
  }
  
  X_control <- cbind(data$propensity, data$X_train)
  X_treat <- cbind(data$propensity, data$X_train)
  treatment <- data$treatment
  true_event_times <- data$true_event_times
  
  n <- length(y_full)
  fold_ids <- sample(rep(1:K, length.out = n))
  
  results_list <- list()
  
  for (i in seq_len(nrow(param_grid))) {
    k_val <- param_grid$k[i]
    local_hp  <- k_val / sqrt(number_of_trees)
    global_hp <- k_val / sqrt(number_of_trees)
    
    cindex_fold <- numeric(K)
    
    for (fold in 1:K) {
      test_idx <- which(fold_ids == fold)
      train_idx <- setdiff(seq_len(n), test_idx)
      
      y_train <- y_full[train_idx]
      status_train <- status_full[train_idx]
      X_train_control <- X_control[train_idx, , drop = FALSE]
      X_train_treat <- X_treat[train_idx, , drop = FALSE]
      treat_train <- treatment[train_idx]
      
      surv_time_test <- exp(y_full[test_idx])
      status_test <- status_full[test_idx]
      X_test_control <- X_control[test_idx, , drop = FALSE]
      X_test_treat <- X_treat[test_idx, , drop = FALSE]
      treat_test <- treatment[test_idx]
      
      fit <- ShrinkageTrees::CausalShrinkageForest(
        y = y_train,
        status = status_train,
        X_train_control = X_train_control,
        X_train_treat = X_train_treat,
        X_test_control = X_test_control,
        X_test_treat = X_test_treat,
        outcome_type = "right-censored",
        treatment_indicator_train = treat_train,
        treatment_indicator_test = treat_test,
        number_of_trees_control = number_of_trees,
        number_of_trees_treat = number_of_trees,
        local_hp_control = local_hp,
        global_hp_control = global_hp,
        local_hp_treat = local_hp,
        global_hp_treat = global_hp,
        prior_type_control = "horseshoe",
        prior_type_treat = "horseshoe",
        N_post = N_post,
        N_burn = N_burn,
        sigma = sigma,
        verbose = FALSE,
        timescale = "log",
        store_posterior_sample = FALSE
      )
      
      pred <- exp(fit$test_predictions)
      
      c_index <- survival::concordance(
        survival::Surv(surv_time_test, status_test) ~ pred
      )$concordance
      
      cindex_fold[fold] <- c_index
    }
    
    results_list[[i]] <- data.frame(
      k = k_val,
      Mean_CIndex = mean(cindex_fold),
      SD_CIndex = sd(cindex_fold)
    )
  }
  
  cv_summary <- do.call(rbind, results_list)
  
  # Select params with highest mean C-index
  best_row <- cv_summary[which.max(cv_summary$Mean_CIndex), ]
  
  local_hp  <- best_row$k / sqrt(number_of_trees)
  global_hp <- best_row$k / sqrt(number_of_trees)
  
  # Refit on full data using best params
  final_fit <- ShrinkageTrees::CausalShrinkageForest(
    y = y_full,
    status = status_full,
    outcome_type = "right-censored",
    X_train_control = X_control,
    X_train_treat = X_treat,
    treatment_indicator_train = treatment,
    number_of_trees_control = number_of_trees,
    number_of_trees_treat = number_of_trees,
    local_hp_control = local_hp,
    global_hp_control = global_hp,
    local_hp_treat = local_hp,
    global_hp_treat = global_hp,
    prior_type_control = "horseshoe",
    prior_type_treat = "horseshoe",
    N_post = N_post,
    N_burn = N_burn,
    sigma = sigma,
    verbose = FALSE,
    timescale = "log",
    store_posterior_sample = TRUE,
  )
  
  # Posterior samples of τ(x)
  cate_samples <- final_fit$train_predictions_sample_treat
  
  # Pointwise quantities
  cate_mean <- final_fit$train_predictions_treat
  cate_rpehe <- sqrt(mean((cate_mean - data$true_cate)^2))
  cate_abs_bias <- mean(abs(cate_mean - data$true_cate))
  
  # Per-observation posterior intervals for τ(x)
  cate_ci <- apply(cate_samples, 2, quantile, probs = c(0.025, 0.975))
  cate_coverage <- mean(cate_ci[1,] <= data$true_cate & cate_ci[2,] >= data$true_cate)
  cate_ci_length <- mean(cate_ci[2,] - cate_ci[1,])
  
  # ATE estimate and uncertainty
  ate_samples <- rowMeans(cate_samples)
  ate_estimate <- mean(ate_samples)
  ate_ci <- quantile(ate_samples, probs = c(0.025, 0.975))
  ate_coverage <- unname(ate_ci[1] <= data$true_ate & ate_ci[2] >= data$true_ate)
  ate_ci_length <- unname(diff(ate_ci))
  
  delta <- 0.25 * sd(data$true_cate)
  cate_detection_power <- compute_cate_detection(
    cate_ci   = cate_ci,
    true_cate = data$true_cate,
    delta     = delta
  )
  
  # RMSE for total outcome prediction
  rmse <- sqrt(mean((final_fit$train_predictions - data$true_event_times)^2))
  
  # Compute the C-index
  if (requireNamespace("survival", quietly = TRUE)) {
    c_index <- survival::concordance(
      survival::Surv(exp(data$follow_up), data$status) ~ exp(final_fit$train_predictions)
    )$concordance
  } else {
    c_index <- NA
  }
  
  return(list(
    postmean_sigma = mean(final_fit$sigma),
    ate = ate_estimate,
    ate_ci = ate_ci,
    ate_coverage = ate_coverage,
    ate_ci_length = ate_ci_length,
    cate_rpehe = cate_rpehe,
    cate_abs_bias  = cate_abs_bias,
    cate_coverage = cate_coverage,
    cate_ci_length = cate_ci_length,
    rmse = rmse,
    cate_detection_power  = cate_detection_power,
    c_index = c_index
  ))
}


run_single_simulation <- function(n_train,
                                  p_feat,
                                  sigma,
                                  s_prog,
                                  s_treat,
                                  cens_scale,
                                  k,
                                  linear,
                                  param_grid,
                                  N_post = 2000,
                                  N_burn = 2000,
                                  n_IV = 0, iv_strength = 0,
                                  n_M = 0, m_strength = 0,
                                  seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  dt <- data_gen(n_train, p_feat, sigma, s_prog, s_treat, cens_scale, 
                 n_IV = n_IV, iv_strength = iv_strength, 
                 n_M = n_M, m_strength = m_strength, # <--- PASS TO GEN
                 linear = linear)
  
  res_CHF <- evaluate_CHF(
    dt,
    number_of_trees_control = 200,
    number_of_trees_treat = 200,
    local_hp_control = k / sqrt(200),
    local_hp_treat = k / sqrt(200),
    global_hp_control = k / sqrt(200),
    global_hp_treat = k / sqrt(200),
    N_post = N_post,
    N_burn = N_burn
  )
  
  res_CHF_CV <- evaluate_CHF_CV(
    dt,
    param_grid = param_grid,
    K = 5,
    number_of_trees = 200,
    N_post = N_post,
    N_burn = N_burn
  )
  
  
  all_methods <- list(
    CHF            = res_CHF,
    CHF_CV         = res_CHF_CV
  )
  
  extract <- function(res) {
    data.frame(
      sigma          = res$postmean_sigma,
      ATE            = res$ate,
      ATE_bias       = res$ate - dt$true_ate,
      ATE_coverage   = res$ate_coverage,
      ATE_CI_Length  = res$ate_ci_length,
      CATE_RPEHE     = res$cate_rpehe,
      CATE_coverage  = res$cate_coverage,
      CATE_CI_Length = res$cate_ci_length,
      CATE_bias = res$cate_abs_bias,
      CATE_Detection_Power  = res$cate_detection_power,
      RMSE           = res$rmse,
      C_Index        = res$c_index
    )
  }
  
  result_df <- do.call(rbind, lapply(all_methods, extract))
  result_df$Method <- names(all_methods)
  rownames(result_df) <- NULL
  result_df <- result_df[, c("Method", setdiff(names(result_df), "Method"))]
  
  return(result_df)
}

run_parallel_simulations <- function(M,
                                     n_train,
                                     p_feat,
                                     sigma,
                                     s_prog,
                                     s_treat,
                                     cens_scale,
                                     k,
                                     linear,
                                     param_grid,
                                     N_post = 2000,
                                     N_burn = 2000) {
  
  results <- foreach(
    i = 1:M,
    .combine = rbind,
    .packages = c("ShrinkageTrees", "survival"),
    .export = c(
      # core components
      "run_single_simulation",
      "data_gen",
      
      # Individual methods
      "evaluate_CHF",
      "evaluate_CHF_CV"
    )
  ) %dopar% {
    
    out <- run_single_simulation(
      n_train = n_train,
      p_feat = p_feat,
      sigma = sigma,
      s_prog = s_prog,
      s_treat = s_treat,
      cens_scale = cens_scale,
      k = k,
      linear = linear,
      param_grid = param_grid,
      N_post = N_post,
      N_burn = N_burn,
      seed = i
    )
    
    cat("Iteration:", i, "done\n")
    
    out
  }
  
  return(results)
}

run_simulations_over_p <- function(
    p_values,
    cens_scales,
    M,
    n_train,
    sigma,
    s_prog,
    s_treat,
    k,
    linear,
    param_grid,
    N_post = 2000,
    N_burn = 2000 
) {
  
  if (length(p_values) != length(cens_scales)) {
    stop("p_values and cens_scales must be the same length.")
  }
  
  all_results <- Map(function(p_val, cs_val) {
    
    cat("Running simulations for p =", p_val, "cens_scale =", cs_val, "\n")
    
    res <- run_parallel_simulations(
      M = M,
      n_train = n_train,
      p_feat = p_val,
      sigma = sigma,
      s_prog = s_prog,
      s_treat = s_treat,
      cens_scale = cs_val,
      k = k,
      linear = linear,
      param_grid = param_grid,
      N_post = N_post,
      N_burn = N_burn
    )
    
    res$p_feat <- p_val
    res$cens_scale <- cs_val
    res
  },
  p_values, cens_scales)
  
  do.call(rbind, all_results)
}




#############################

data("pdac")

time      <- pdac$time
status    <- pdac$status
treatment <- pdac$treatment
covariates <- pdac[, !(colnames(pdac) %in% c("time", "status", "treatment"))]
covariates <- scale(covariates)
n_train <- nrow(covariates)
p_feat  <- ncol(covariates)



data_gen <- function(n_train, p_feat, sigma, s_prog, s_treat, cens_scale, 
                     n_IV = 0, iv_strength = 0, # Existing IV args
                     n_M = 0, m_strength = 0,   # <--- NEW M-BIAS ARGUMENTS
                     linear = NULL) {
  
  # Start with original covariates
  X_train_original <- covariates
  treatment_i <- treatment
  n <- nrow(X_train_original)
  p <- ncol(X_train_original)
  
  # Initialize the combined matrix
  X_generated <- NULL
  
  # --- 1. GENERATE Z-VARIABLES (IVs) ---
  if (n_IV > 0) {
    trt_centered <- treatment_i - 0.5 
    Z_noise <- matrix(rnorm(n * n_IV), nrow = n, ncol = n_IV)
    Z_signal <- matrix(trt_centered * iv_strength, nrow = n, ncol = n_IV)
    Z_matrix <- scale(Z_noise + Z_signal)
    colnames(Z_matrix) <- paste0("Z_IV", 1:n_IV)
    X_generated <- cbind(X_generated, Z_matrix)
  }
  
  # --- 2. GENERATE M-VARIABLES (Colliders) ---
  # Structure: A <--- U1 ---> M <--- U2 ---> Y
  
  U2_effect_on_Y <- rep(0, n) # Default 0
  
  if (n_M > 0) {
    # U1: Unmeasured confounder causing A (correlated with fixed treatment)
    trt_centered <- treatment_i - 0.5
    # We simulate U1 as being "caused" by A here to force the correlation
    U1 <- trt_centered + rnorm(n, 0, 1) 
    
    # U2: Unmeasured confounder causing Y
    U2 <- matrix(rnorm(n * n_M), nrow = n, ncol = n_M)
    
    # M: The Collider (Observed) = U1 + U2
    # We sum U1 (vector) with every column of U2
    M_matrix <- U2 + as.vector(U1)
    
    # Scale M for the model
    M_matrix <- scale(M_matrix)
    colnames(M_matrix) <- paste0("M_Collider", 1:n_M)
    
    # Append M to generated features
    if (is.null(X_generated)) {
      X_generated <- M_matrix
    } else {
      X_generated <- cbind(X_generated, M_matrix)
    }
    
    # Calculate the effect of U2 on Y (This is the "Bias Path")
    # If m_strength > 0, U2 drives Y. 
    # Since M is in the model, the model will use M to proxy U2.
    # But M is linked to U1 (and A), so this creates bias.
    # We sum the effects of all U2 columns
    U2_effect_on_Y <- rowSums(U2) * m_strength
  }
  
  # Combine Original X with Generated Variables (Z or M)
  if (!is.null(X_generated)) {
    X_train_combined <- cbind(X_train_original, X_generated)
  } else {
    X_train_combined <- X_train_original
  }
  
  # --- 3. GENERATE OUTCOMES ---
  
  beta_f <- rnorm(p, 0, 2) * rbinom(p, 1, s_prog)
  f_x <- 6 + as.numeric(beta_f %*% t(X_train_original))
  
  beta_tau <- rnorm(p, 1, 0.5) * rbinom(p, 1, s_treat)
  tau_x <- 2 + as.numeric(beta_tau %*% t(X_train_original))
  
  true_ate_raw <- sum(beta_tau * colMeans(X_train_original))
  
  # Log-normal survival times
  # NOTE: We add U2_effect_on_Y here!
  log_T <- f_x + treatment_i * tau_x + U2_effect_on_Y + rnorm(n, 0, sigma)
  
  T <- exp(log_T)
  C <- rexp(n, rate = cens_scale)
  Y <- pmin(T, C)
  status <- as.numeric(T <= C)
  sd_un <- sd(log_T)
  
  
  prop_fit <- ShrinkageTrees::HorseTrees(
    y = treatment_i,
    X_train = X_train_combined,
    outcome_type = "binary",
    k = 0.1,
    N_post = 3000,
    N_burn = 2000,
    verbose = FALSE
  )
  propensity_i <- pnorm(prop_fit$train_predictions)
  
  
  return(list(
    X_train                = X_train_combined, 
    treatment              = treatment_i,
    propensity             = propensity_i,
    follow_up              = log(Y) / sd_un,
    status                 = status,
    true_event_times       = log(T) / sd_un,
    uncensored_event_times = log(T) / sd_un,
    true_cate              = tau_x / sd_un,
    true_ate               = true_ate_raw / sd_un,
    sample_ate             = mean(tau_x) / sd_un,
    obs_sigma              = sigma / sd_un
  ))
}

# ... (Keep all your function definitions: compute_cate_detection, evaluate_CHF, evaluate_CHF_CV, run_single_simulation, data_gen) ...

# ==============================================================================
# EXECUTION START
# ==============================================================================

# Retrieve command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  num_cores <- as.integer(args[1]) - 1
} else {
  num_cores <- parallel::detectCores() - 1
}

registerDoParallel(cores = num_cores)

cat("Number of cores being used (1 free):", num_cores, "\n")
cat("SIMULATION: PDAC Sensitivity Analysis (IV & M-Bias)\n")

# --- Simulation Parameters ---
M <- 3*num_cores
n_train <- nrow(covariates)
N_post <- 3000
N_burn <- 2000
s_prog <- 0.1
s_treat <- 0.05
sigma <- 3/4
cens_scale_hat <- 1/400
p_vals <- c(ncol(covariates))

# Define Hyperparameter Grid for CV
param_grid_for_CV <- expand.grid(
  k = c(0.05, 0.1, 0.25, 0.5, 0.75, 1)
)

scenarios <- list(
  # 1. Baseline (Reference)
  list(name = "01_Baseline",         n_IV=0, iv_strength=0, n_M=0,  m_strength=0),
  
  # --- TEST A: THE "BAIT" TEST (Variable Strength, Fixed Low Dim) ---
  # n_M is small (5), so the "minefield" is sparse.
  # We test how "attractive" M must be (as a proxy for U2) to trick the model.
  list(name = "02_M_Weak",           n_IV=0, iv_strength=0, n_M=5,  m_strength=0.5),
  list(name = "03_M_Moderate",       n_IV=0, iv_strength=0, n_M=5,  m_strength=1.0),
  list(name = "04_M_Strong",         n_IV=0, iv_strength=0, n_M=5,  m_strength=3.0),
  list(name = "05_M_Extreme",        n_IV=0, iv_strength=0, n_M=5,  m_strength=5.0),
  
  # --- TEST B: THE "MINEFIELD" TEST (Variable Dim, Fixed Moderate Strength) ---
  # Strength is fixed at 1.0 (moderate bait).
  # We test if the model fails simply because there are *so many* colliders to accidentally pick.
  # Note: Since your code sums U2 effects, higher n_M also implicitly adds more variance to Y.
  list(name = "06_Dim_Medium",       n_IV=0, iv_strength=0, n_M=10, m_strength=1.0),
  list(name = "07_Dim_High",         n_IV=0, iv_strength=0, n_M=25, m_strength=1.0),
  list(name = "08_Dim_Massive",      n_IV=0, iv_strength=0, n_M=50, m_strength=1.0),
  
  # --- TEST C: THE "TORTURE" TEST (High Dim + High Strength) ---
  # Many variables, and all of them are highly predictive of Y (strong bait).
  # If the forest selects these M variables, it will induce significant bias.
  list(name = "09_Stress_Med",       n_IV=0, iv_strength=0, n_M=10, m_strength=3.0),
  list(name = "10_Stress_High",      n_IV=0, iv_strength=0, n_M=25, m_strength=3.0),
  list(name = "11_Stress_Max",       n_IV=0, iv_strength=0, n_M=50, m_strength=3.0)
)

all_scenario_results <- list()

# --- Main Loop ---
for (scen in scenarios) {
  cat(sprintf("\n--- Running Scenario: %s ---\n", scen$name))
  
  # Export necessary objects to workers
  # Note: param_grid_for_CV and cens_scale_hat added to export
  res <- foreach(i = 1:M, .combine = rbind, 
                 .packages = c("ShrinkageTrees", "survival"), 
                 .export = c("run_single_simulation", "data_gen", 
                             "evaluate_CHF", "evaluate_CHF_CV", "compute_cate_detection",
                             "covariates", "treatment", 
                             "param_grid_for_CV", "cens_scale_hat")) %dopar% {
                               
                               run_single_simulation(
                                 n_train = n_train, 
                                 p_feat = p_feat, 
                                 sigma = sigma,
                                 s_prog = s_prog, 
                                 s_treat = s_treat, 
                                 cens_scale = cens_scale_hat,
                                 k = 0.3, 
                                 linear = NULL, 
                                 param_grid = param_grid_for_CV,
                                 N_post = N_post, 
                                 N_burn = N_burn,
                                 
                                 # Sensitivity Args
                                 n_IV = scen$n_IV, 
                                 iv_strength = scen$iv_strength,
                                 n_M = scen$n_M,           
                                 m_strength = scen$m_strength, 
                                 
                                 seed = i
                               )
                             }
  
  res$Scenario <- scen$name
  all_scenario_results[[scen$name]] <- res
}

# --- Combine and Save ---

# Combine list of dataframes into one large dataframe
final_combined_df <- do.call(rbind, all_scenario_results)

# Define output file path 
output_file <- file.path(Sys.getenv("TMPDIR"), "PDAC_sensitivity_M_prop_output.rds")

# Print message
cat("Saving all settings results to:", output_file, "\n")

# Save the CORRECT object
saveRDS(final_combined_df, file = output_file)

# Confirm successful save
cat("All results successfully saved.\n")






# results <- readRDS("/Users/tijnjacobs/Library/CloudStorage/OneDrive-VrijeUniversiteitAmsterdam/Documents/GitHub/ShrinkageTrees/simulations revision/Data analysis simulations/PDAC_sensitivity_M_prop_output.rds")
# 
# summary_table <- results %>%
#   group_by(Method, Scenario) %>% # <--- CRITICAL ADDITION
#   summarise(
#     # ... (keep your existing summary stats) ...
#     CATE_abs_bias  = mean(abs(CATE_bias), na.rm = TRUE),
#     CATE_RPEHE       = mean(CATE_RPEHE, na.rm = TRUE),
#     CATE_coverage   = mean(CATE_coverage, na.rm = TRUE),
#     # ...
#     .groups = "drop"
#   ) %>%
#   arrange(Method, match(Scenario, c("Baseline", "Z_Bias_Moderate", "Z_Bias_Strong"))) # Optional: Sort logically
# 
# print(summary_table)
