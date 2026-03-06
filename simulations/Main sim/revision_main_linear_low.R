### Simulation deeper Friedman ###

library(ShrinkageTrees)
library(foreach)
library(doParallel)

source("/home/tjacobs/evaluation_functions_main.R")

data_gen <- function(n_train, p_feat, sigma, s_prog, s_treat, cens_scale = 1, linear) {
  X_train <- matrix(runif(n_train * p_feat), n_train, p_feat)
  
  propensity <- pnorm(-0.5 + 0.4 * X_train[,1] - 0.1 * X_train[,3] + 0.3 * X_train[,5])
  treatment <- rbinom(n_train, 1, propensity)
  
  beta_prog <- rnorm(p_feat, 0, 1) * rbinom(p_feat, 1, s_prog)
  mu <- beta_prog %*% t(X_train)
  
  if (linear) {
    
    beta_treat <- rnorm(p_feat, 0, 1) * rbinom(p_feat, 1, s_treat)
    
    tau <- 1 + X_train[, 1] - 2 * X_train[, 2] + 3 * X_train[, 3] - 
      4 * X_train[, 4] + 5 * X_train[, 5] + beta_treat %*% t(X_train)
    
    true_ate <- 1 - 2 * 1/2 + 3 * 1/2 - 4 * 1/2 + 5 * 1/2 + sum(beta_prog) * 1/2
    
  } else {
    
    tau <- 10 * sin(pi * X_train[, 1] * X_train[, 2]) +
      20 * (X_train[, 3] - 0.5)^2 +
      10 * X_train[, 4] +
      5  * X_train[, 5]
    
    true_ate <- 10 * 0.5246745 + 20 * 1/12 + 10 * 1/2 + 5 * 1/2
  }
  
  true_event_times <- mu + (treatment - 0.5) * tau
  uncensored_event_times <- true_event_times + rnorm(n_train, 0, sigma)
  sd_un <- sd(uncensored_event_times)
  uncensored_event_times <- uncensored_event_times / sd_un
  true_event_times <- true_event_times / sd_un
  
  C <- log(rexp(n_train, cens_scale)) + min(uncensored_event_times)
  follow_up <- pmin(uncensored_event_times, C)
  status <- as.numeric(uncensored_event_times <= C)
  
  return(list(
    X_train = X_train,
    treatment = treatment,
    propensity = propensity,
    follow_up = as.numeric(follow_up),
    status = status,
    true_event_times = as.numeric(true_event_times),
    uncensored_event_times = as.numeric(uncensored_event_times),
    true_cate = as.vector(tau) / sd_un,
    true_ate = true_ate / sd_un,
    sample_ate = mean(tau) / sd_un,
    obs_sigma = sigma / sd_un
  ))
}


find_cens_scale <- function(p,
                            n_train = 100,
                            sigma = sqrt(3),
                            s_prog,
                            s_treat,
                            target_cens,
                            linear,
                            M = 10000,
                            verbose = FALSE) {
  
  # Inner function: return absolute difference between actual and target censoring
  estimate_cens_diff <- function(cens_scale) {
    rates <- replicate(M, {
      dt <- data_gen(n_train, p, sigma, s_prog, s_treat, cens_scale, linear)
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
cat("SIMULATION: revision main linear\n")

M <- 1000
n_train <- 200
sigma <- sqrt(3)
s_prog <- 0.1
s_treat <- 0.05
N_post <- 5000
N_burn <- 2500

param_grid_for_CV <- expand.grid(
  k = c(0.25, 0.3, 0.35)
)

p_vals <- c(100, 1000, 5000)
# cens_scales_linear <- sapply(p_vals, function(p) find_cens_scale(p = p, n_train = n_train, sigma = sigma, s_prog = s_prog, s_treat = s_treat, target_cens = 0.35, linear = TRUE, M = 500))
cens_scales_linear <- c(0.02470207, 0.02403064, 0.02390695)

results_linear <- run_simulations_over_p(
  p_values = p_vals,
  cens_scales = cens_scales_linear,
  M = M,
  n_train = n_train,
  sigma = sigma,
  s_prog = s_prog,
  s_treat = s_treat,
  k = 0.1, 
  linear = TRUE,
  param_grid = param_grid_for_CV,
  N_post = N_post,
  N_burn = N_burn
)

info_object <- list(
  simulation_name = "revision_main_linear_friedman",
  date = Sys.time(),
  
  data_generating_process = list(
    n_train   = n_train,
    p_values  = p_vals,
    sigma     = sigma,
    s_prog    = s_prog,
    s_treat   = s_treat,
    target_censoring = 0.35,
    cens_scales = cens_scales_linear
  ),
  
  model_settings = list(
    N_post = N_post,
    N_burn = N_burn,
    k_fixed = 0.1,
    k_grid  = param_grid_for_CV$k
  ),
  
  simulation_settings = list(
    M = M,
    num_cores = num_cores,
    parallel_backend = "doParallel"
  ),
  
  software = list(
    R_version = R.version.string,
    package = "ShrinkageTrees"
  ),
  
  notes = "Linear Friedman-style DGP with high-dimensional covariates; censoring calibrated per p."
)


combined_results <- list(info_object = info_object,
                         results_linear = results_linear)

# Define output file path 
# (! NAME MUST BE FILENAME_output.rds !)
output_file <- file.path(Sys.getenv("TMPDIR"), "revision_main_linear_low_output.rds")

# Print message
cat("Saving all settings results to:", output_file, "\n")

# Save the combined list
saveRDS(combined_results, file = output_file)

# Confirm successful save
cat("All results successfully saved in one file.\n")









# 
# 
# library(dplyr)
# 
# # Load results
# res <- readRDS(
#   "~/Library/CloudStorage/OneDrive-VrijeUniversiteitAmsterdam/Documents/GitHub/ShrinkageTrees/revision/Main sim/revision_main_linear_low_output.rds"
# )
# 
# # Extract and summarize CATE metrics
# cate_summary <- res$results_linear %>%
#   select(
#     Method, p_feat,
#     CATE_RPEHE,
#     CATE_AbsBias,
#     CATE_coverage,
#     CATE_CI_Length,
#     CATE_Detection_Power0.25,
#     CATE_Detection_Power0.1,
#     CATE_winkler,
#     CATE_csd,
#     CATE_empirical_power
#   ) %>%
#   group_by(Method, p_feat) %>%
#   summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)),
#             .groups = "drop")
# 
# print(cate_summary, n=100)
# 
# 
