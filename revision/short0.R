library(ShrinkageTrees)
library(foreach)
library(doParallel)

source("/home/tjacobs/evaluation_functions_new.R")

data_gen <- function(n_train, p_feat, sigma, s_prog, s_treat, cens_scale, linear = NULL) {
  X_train <- matrix(runif(n_train * p_feat), n_train, p_feat)
  
  # Propensity & treatment
  propensity <- pnorm(-0.5 + 0.4 * X_train[, 1])
  treatment  <- rbinom(n_train, 1, propensity)
  
  # Prognostic component
  beta_prog <- rnorm(p_feat, 0, 1) * sort(rbinom(p_feat, 1, s_prog), TRUE)
  mu <- as.numeric(beta_prog %*% t(X_train))
  
  # Baseline nonlinear tau (Friedman-type)
  tau <- 10 * sin(pi * X_train[, 1] * X_train[, 2]) +
    20 * (X_train[, 3] - 0.5)^2 +
    10 * X_train[, 4] +
    5  * X_train[, 5]
  
  # Analytic ATE for the baseline part (same as your original)
  base_ate <- 10 * 0.5246745 + 20 * 1/12 + 10 * 1/2 + 5 * 1/2
  
  ## Add sparse extra main & interaction terms to tau
  m_extra <- max(1, round(p_feat * s_treat))  # total extra terms
  
  # Choose split between main and interaction terms (can adjust)
  m_main <- min(p_feat, floor(m_extra / 2))
  m_int2 <- m_extra - m_main
  
  extra_ate <- 0  # analytic contribution to E[tau(X)] from extra terms
  
  # Main effects
  if (m_main > 0) {
    idx_main   <- sample.int(p_feat, m_main)
    beta_main  <- rnorm(m_main, mean = 0, sd = 1)
    
    # Add to tau
    tau <- tau + as.numeric(X_train[, idx_main, drop = FALSE] %*% beta_main)
    
    # For X_j ~ Uniform(0,1), E[X_j] = 1/2
    extra_ate <- extra_ate + sum(beta_main) * 1/2
  }
  
  # Pairwise interactions
  if (m_int2 > 0) {
    # Draw m_int2 random pairs (j, k)
    pair_idx <- replicate(m_int2, sample.int(p_feat, 2), simplify = TRUE)  # 2 x m_int2
    beta_int <- rnorm(m_int2, mean = 0, sd = 1)
    
    for (k in seq_len(m_int2)) {
      j1 <- pair_idx[1, k]
      j2 <- pair_idx[2, k]
      tau <- tau + beta_int[k] * (X_train[, j1] * X_train[, j2])
    }
    
    # For independent X_j, X_k ~ U(0,1): E[X_j X_k] = 1/4
    extra_ate <- extra_ate + sum(beta_int) * 1/4
  }
  
  # Total analytic ATE before standardisation
  true_ate_raw <- base_ate + extra_ate
  
  ## Generate event times and censoring
  
  true_event_times <- mu + (treatment - 0.5) * tau
  uncensored_event_times <- true_event_times + rnorm(n_train, 0, sigma)
  
  sd_un <- sd(uncensored_event_times)
  uncensored_event_times <- uncensored_event_times / sd_un
  true_event_times       <- true_event_times / sd_un
  
  # Exponential censoring on log scale
  C <- log(rexp(n_train, cens_scale)) + min(uncensored_event_times)
  follow_up <- pmin(uncensored_event_times, C)
  status    <- as.numeric(uncensored_event_times <= C)
  
  return(list(
    X_train               = X_train,
    treatment             = treatment,
    propensity            = propensity,
    follow_up             = as.numeric(follow_up),
    status                = status,
    true_event_times      = as.numeric(true_event_times),
    uncensored_event_times = as.numeric(uncensored_event_times),
    true_cate             = as.vector(tau) / sd_un,
    true_ate              = true_ate_raw / sd_un,
    sample_ate            = mean(tau) / sd_un,
    obs_sigma             = sigma / sd_un
  ))
}


find_cens_scale <- function(p,
                            n_train,
                            sigma,
                            s_prog,
                            s_treat,
                            target_cens,
                            M = 1000,
                            verbose = FALSE) {
  
  # Inner function: return absolute difference between actual and target censoring
  estimate_cens_diff <- function(cens_scale) {
    rates <- replicate(M, {
      dt <- data_gen(n_train, p, sigma, s_prog, s_treat, cens_scale, linear = NULL)
      mean(1 - dt$status)
    })
    diff <- abs(mean(rates) - target_cens)
    if (verbose) cat("cens_scale =", round(cens_scale, 4), 
                     "| estimated cens =", round(mean(rates), 3), "\n")
    return(diff)
  }
  
  # Optimize over a reasonable range
  opt_result <- optimize(estimate_cens_diff, interval = c(0.001, 2), tol = 1e-3)
  
  return(opt_result$minimum)
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

M <- num_cores
n_train <- 100
sigma <- 1
s_prog <- 0.1
s_treat <- 0.05
N_post <- 3
N_burn <- 2

p_vals <- c(100, 1000, 10000)
cens_scales_medium <- c(0.05453025, 0.03694038, 0.02376848)


param_grid_for_CV <- expand.grid(
  k = c(0.05, 0.1, 0.25, 0.5, 0.75, 1)
)


cat("\nStarting with medium censoring scenario\n\n")
results_medium <- run_simulations_over_p(
  p_values = p_vals,
  cens_scales = cens_scales_medium,
  M = M,
  n_train = n_train,
  sigma = sigma,
  s_prog = s_prog,
  s_treat = s_treat,
  param_grid = param_grid_for_CV,
  k = 0.1,
  N_post = N_post,
  N_burn = N_burn
)

combined_results <- list(results_medium = results_medium)


# Define output file path 
# (! NAME MUST BE FILENAME_output.rds !)
output_file <- file.path(Sys.getenv("TMPDIR"), "short0_output.rds")

# Print message
cat("Saving all settings results to:", output_file, "\n")

# Save the combined list
saveRDS(combined_results, file = output_file)

# Confirm successful save
cat("All results successfully saved in one file.\n")






library(dplyr)
library(tidyr)

res <- readRDS(
  "/Users/tijnjacobs/Library/CloudStorage/OneDrive-VrijeUniversiteitAmsterdam/Documents/GitHub/ShrinkageTrees/simulations revision/short0_output.rds"
)$results_medium

rename_methods <- function(x) {
  recode(x,
         "CHF (k=0.1)" = "CHF",
         "DART_AFT"    = "AFT_DART",
         "S_BCF_AFT"   = "AFT_S_BCF",
         "IndivAFT"    = "AFT_BART",
         "BCF_AFT"     = "AFT_BCF",
         "CHF_CV"      = "CHF_CV")
}

res <- res %>%
  mutate(Method = rename_methods(Method))


cate_table <- res %>%
  group_by(Method, p_feat) %>%
  summarise(
    # Accuracy
    CATE_RMSE              = mean(CATE_RPEHE, na.rm = TRUE),
    CATE_AbsBias           = mean(CATE_AbsBias, na.rm = TRUE),

    # Uncertainty
    CATE_Coverage          = mean(CATE_coverage, na.rm = TRUE),
    CATE_CI_Length         = mean(CATE_CI_Length, na.rm = TRUE),

    # Significance-based metrics
    Empirical_Power        = mean(CATE_empirical_power, na.rm = TRUE),
    CSD                    = mean(CATE_csd, na.rm = TRUE),
    Winkler_Score          = mean(CATE_winkler, na.rm = TRUE),

    # Thresholded detection power
    DetectionPower_0.10    = mean(CATE_Detection_Power0.1, na.rm = TRUE),
    DetectionPower_0.25    = mean(CATE_Detection_Power0.25, na.rm = TRUE),

    .groups = "drop"
  ) %>%
  arrange(p_feat, Method)

cate_table
