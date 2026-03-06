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
                                  n_IV = 0,       
                                  iv_strength = 0,
                                  seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  dt <- data_gen(n_train, p_feat, sigma, s_prog, s_treat, cens_scale, 
                 n_IV = n_IV, iv_strength = iv_strength, linear = linear)
  
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
  
  # res_CHF_CV <- evaluate_CHF_CV(
  #   dt,
  #   param_grid = param_grid,
  #   K = 5,
  #   number_of_trees = 200,
  #   N_post = N_post,
  #   N_burn = N_burn
  # )
  
  
  all_methods <- list(
    # CHF_CV         = res_CHF_CV,
    CHF            = res_CHF
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
                     n_IV = 0, iv_strength = 0, # <--- NEW ARGUMENTS
                     linear = NULL) {
  
  # Start with original covariates
  X_train_original <- covariates
  treatment_i <- treatment
  
  n <- nrow(X_train_original)
  p <- ncol(X_train_original)
  
  # --- 1. GENERATE Z-VARIABLES (IVs) ---
  if (n_IV > 0) {
    # Create Z based on Treatment (Reverse Causality for simulation)
    # Z = Noise + (Strength * Treatment)
    
    # We center treatment (0/1 -> -0.5/0.5) to keep Z roughly centered around 0
    trt_centered <- treatment_i - 0.5 
    
    # Generate random normal noise for Z
    Z_noise <- matrix(rnorm(n * n_IV), nrow = n, ncol = n_IV)
    
    # Add the treatment signal to Z
    # Stronger 'iv_strength' = Higher correlation between Z and A
    Z_signal <- matrix(trt_centered * iv_strength, nrow = n, ncol = n_IV)
    
    Z_matrix <- Z_noise + Z_signal
    
    # Optional: Scale Z so it looks like other standardized covariates
    Z_matrix <- scale(Z_matrix)
    colnames(Z_matrix) <- paste0("Z_IV", 1:n_IV)
    
    # Combine original X with new Z
    X_train_combined <- cbind(X_train_original, Z_matrix)
    
  } else {
    X_train_combined <- X_train_original
  }
  
  # --- 2. GENERATE OUTCOMES (Using ONLY X, NOT Z) ---
  
  # Prognostic function f(x) - Uses ORIGINAL covariates only
  # We do not want Z to affect the outcome (that would be confounding, not Z-bias)
  beta_f <- rnorm(p, 0, 2) * rbinom(p, 1, s_prog)
  f_x <- 6 + as.numeric(beta_f %*% t(X_train_original))
  
  # Treatment effect tau(x) - Uses ORIGINAL covariates only
  beta_tau <- rnorm(p, 1, 0.5) * rbinom(p, 1, s_treat)
  tau_x <- 2 + as.numeric(beta_tau %*% t(X_train_original))
  
  true_ate_raw <- sum(beta_tau * colMeans(X_train_original))
  
  # Log-normal survival times
  log_T <- f_x + treatment_i * tau_x + rnorm(n, 0, sigma)
  T <- exp(log_T)
  
  # Exponential censoring
  C <- rexp(n, rate = cens_scale)
  
  Y <- pmin(T, C)
  status <- as.numeric(T <= C)
  
  sd_un <- sd(log_T)
  
  
  
  # Refit propensity on augmented X
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
    X_train                = X_train_combined, # This now includes Z!
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


# Retrieve command-line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) > 0) {
  num_cores <- as.integer(args[1]) - 1
} else {
  num_cores <- parallel::detectCores() - 1
}

registerDoParallel(cores = num_cores)

cat("Number of cores being used (1 free):", num_cores, "\n")
cat("SIMULATION: DEEPER FRIEDMAN2 high\n")

M <- 1000
n_train <- nrow(covariates)
N_post <- 3000
N_burn <- 2000
s_prog <- 0.1
s_treat <- 0.05
sigma <- 3/4
cens_scale_hat <- 1/400
p_vals <- c(ncol(covariates))
cens_scales_medium <- c(cens_scale_hat)

param_grid_for_CV <- expand.grid(
  k = c(0.05, 0.1, 0.25, 0.5, 0.75, 1)
)

scenarios <- list(
  # 1. Baseline (Reference)
  list(name = "01_Baseline",         n_IV = 0,  iv_strength = 0),
  
  # --- TEST A: INCREASING STRENGTH (Fixed Low Dimension n=5) ---
  # How "instrumental" does the variable have to be to confuse the model?
  # list(name = "02_Str_Weak",         n_IV = 5,  iv_strength = 0.5), # ~20% correlation
  # list(name = "03_Str_Moderate",     n_IV = 5,  iv_strength = 1.0), # ~45% correlation
  # list(name = "04_Str_High",         n_IV = 5,  iv_strength = 3.0), # ~80% correlation
  # list(name = "05_Str_Extreme",      n_IV = 5,  iv_strength = 5.0), # ~90%+ correlation (Almost collinear)
  
  # --- TEST B: INCREASING DIMENSION (Fixed Moderate Strength = 1.0) ---
  # Can the sparsity prior handle a flood of moderately correlated noise?
  list(name = "06_Dim_Medium",       n_IV = 10, iv_strength = 1.0),
  # list(name = "07_Dim_High",         n_IV = 25, iv_strength = 1.0),
  # list(name = "08_Dim_VeryHigh",     n_IV = 50, iv_strength = 1.0),
  
  # --- TEST C: THE "STRESS TEST" (High Dim + High Strength) ---
  # The hardest settings: Many variables that all look like the treatment.
  # list(name = "09_Stress_Med",       n_IV = 10, iv_strength = 3.0),
  # list(name = "10_Stress_High",      n_IV = 25, iv_strength = 3.0),
  list(name = "11_Stress_Massive",   n_IV = 50, iv_strength = 3.0)
)

# Initialize storage
all_scenario_results <- list()

for (scen in scenarios) {
  
  cat(sprintf("\n--- Running Scenario: %s (n_IV=%d, strength=%.1f) ---\n", 
              scen$name, scen$n_IV, scen$iv_strength))
  
  # Need to update run_parallel_simulations to accept n_IV/strength 
  # OR just call run_single inside a loop here if you prefer. 
  # Assuming you update run_parallel_simulations similar to run_single_simulation:
  
  res <- foreach(i = 1:M, .combine = rbind, 
                 .packages = c("ShrinkageTrees", "survival"), 
                 .export = c("run_single_simulation", "data_gen", "evaluate_CHF", "evaluate_CHF_CV", "covariates", "treatment")) %dopar% {
                   
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
                     n_IV = scen$n_IV,          # <--- PASS SCENARIO ARGS
                     iv_strength = scen$iv_strength,
                     seed = i
                   )
                 }
  
  # Add scenario label column
  res$Scenario <- scen$name
  res$n_IV <- scen$n_IV
  res$IV_Strength <- scen$iv_strength
  
  all_scenario_results[[scen$name]] <- res
}

# --- Combine all scenarios into one data frame ---
final_combined_df <- do.call(rbind, all_scenario_results)

# --- Define Output File Path ---
# (! NAME MUST BE FILENAME_output.rds !)
output_file <- file.path(Sys.getenv("TMPDIR"), "PDAC_sensitivity_Z_prop_output.rds")

cat("Saving all settings results to:", output_file, "\n")

# --- Save the Correct Object ---
saveRDS(final_combined_df, file = output_file)

cat("All results successfully saved in one file.\n")



# 
# 
# 
# results <- readRDS("/Users/tijnjacobs/Library/CloudStorage/OneDrive-VrijeUniversiteitAmsterdam/Documents/GitHub/ShrinkageTrees/simulations revision/Sensitivity/PDAC_sensitivity_Z_prop_output.rds")
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
