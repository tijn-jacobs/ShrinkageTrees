### Some simulation to practice ###

library(foreach)
library(doParallel)
library(ShrinkageTrees)

source("evaluation_functions.R")


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
                            n_train = 200,
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
cat("SIMULATION: main both\n")

M <- 1
n_train <- 200
sigma <- sqrt(3)
s_prog <- 0.1
s_treat <- 0.05
N_post <- 5
N_burn <- 2

param_grid_for_CV <- expand.grid(
  k = c(0.05, 0.1, 0.25, 0.5, 0.75, 1)
)

p_vals <- c(100, 1000, 5000)

# cens_scales_nonlinear <- sapply(p_vals, function(p) find_cens_scale(p, target_cens = 0.35, s_prog = s_prog, s_treat = s_treat, M = 2000, linear = FALSE))
cens_scales_nonlinear <- c(0.07026994, 0.05060948, 0.03609769)

cat("\nStarting with nonlinear scenario\n\n")
results_nonlinear <- run_simulations_over_p(
  p_values = p_vals,
  cens_scales = cens_scales_nonlinear,
  M = M,
  n_train = n_train,
  sigma = sigma,
  s_prog = s_prog,
  s_treat = s_treat,
  k = 0.1, 
  linear = FALSE,
  param_grid = param_grid_for_CV,
  N_post = N_post,
  N_burn = N_burn
)

combined_results <- list(results_nonlinear = results_nonlinear)


# Define output file path 
# (! NAME MUST BE FILENAME_output.rds !)
output_file <- file.path(Sys.getenv("TMPDIR"), "main_nonlinear_output.rds")

# Print message
cat("Saving all settings results to:", output_file, "\n")

# Save the combined list
saveRDS(combined_results, file = output_file)

# Confirm successful save
cat("All results successfully saved in one file.\n")




