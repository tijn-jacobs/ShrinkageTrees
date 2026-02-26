library(ShrinkageTrees)
library(foreach)
library(doParallel)

source("/home/tjacobs/evaluation_functions_ric.R")



data_gen <- function(n_train,
                     p_feat,
                     par,
                     sigma = 1,
                     eta,
                     seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  qs = par$qs
  qw = par$qw
  cWe = par$cWe
  cWf = par$cWf
  cSe = par$cSe
  cSf = par$cSf
  tau0 = par$tau0
  hetero = par$hetero
  
  ## Covariates
  X_train <- matrix(runif(n_train * p_feat), n_train, p_feat) - 1/2
  
  ## Index sets
  S <- seq_len(qs)
  W <- seq(qs + 1, qs + qw)
  N <- seq(qs + qw + 1, p_feat)
  
  ## Propensity score & treatment
  k <- qs + qw
  scaling_factor <- 1 / sqrt(k)
  
  alpha0 <- 0
  alpha_S <- rep(cSe, qs)
  alpha_W <- rep(cWe, qw)
  
  linpred_e <- alpha0 +
    scaling_factor * (
      X_train[, S, drop = FALSE] %*% alpha_S +
        X_train[, W, drop = FALSE] %*% alpha_W
    )
  
  propensity <- pnorm(linpred_e)
  treatment  <- rbinom(n_train, 1, propensity)
  hist(propensity)
  
  ## Prognostic component
  beta_S <- rep(cSf, qs)
  beta_W <- rep(cWf, qw)
  
  mu <- as.numeric(
    X_train[, S, drop = FALSE] %*% beta_S +
      X_train[, W, drop = FALSE] %*% beta_W
  )
  
  ## Treatment effect
  if (!hetero) {
    tau <- rep(tau0, n_train)
  } else {
    tau_S <- rep(0.5, qs)
    tau <- as.numeric(tau0 + X_train[, S, drop = FALSE] %*% tau_S)
  }
  
  ## Generate log event times
  true_event_times <- mu + (treatment - 0.5) * tau
  uncensored_event_times <- true_event_times +
    rnorm(n_train, 0, sqrt(sigma))
  
  ## Rescale (as in original data_gen)
  sd_un <- sd(uncensored_event_times)
  uncensored_event_times <- uncensored_event_times / sd_un
  true_event_times       <- true_event_times / sd_un
  
  ## Exponential censoring on log scale
  C <- log(rexp(n_train, rate = eta)) + min(uncensored_event_times)
  
  follow_up <- pmin(uncensored_event_times, C)
  status    <- as.numeric(uncensored_event_times <= C)
  
  return(list(
    X_train                = X_train,
    treatment              = treatment,
    propensity             = propensity,
    follow_up              = as.numeric(follow_up),
    status                 = status,
    true_event_times       = as.numeric(true_event_times),
    uncensored_event_times = as.numeric(uncensored_event_times),
    true_cate              = as.vector(tau) / sd_un,
    obs_sigma              = sqrt(sigma) / sd_un,
    sets                   = list(S = S, W = W, N = N),
    params                 = list(
      n_train = n_train,
      p_feat  = p_feat,
      qs      = qs,
      qw      = qw,
      cWe      = cWe,
      cWf      = cWf,
      tau0   = tau0,
      hetero = hetero,
      eta    = eta
    )
  ))
}

compute_cens_scales_parallel <- function(qw_values,
                                         n_train,
                                         p_feat,
                                         par,
                                         sigma,
                                         target_cens = 0.35,
                                         M = 500,
                                         eta_range = c(1e-4, 2),
                                         verbose = TRUE) {
  
  foreach(qw_val = qw_values,
          .combine = c,
          .packages = c("ShrinkageTrees")) %dopar% {
            
            par_loc <- par
            par_loc$qw <- qw_val
            
            if (verbose) {
              cat("Calibrating eta for qw =", qw_val, "\n")
            }
            
            objective <- function(eta) {
              cens_rate <- mean(replicate(M, {
                dt <- data_gen(
                  n_train = n_train,
                  p_feat  = p_feat,
                  par     = par_loc,
                  sigma   = sigma,
                  eta     = eta
                )
                mean(1 - dt$status)
              }))
              
              abs(cens_rate - target_cens)
            }
            
            optimize(objective,
                     interval = eta_range,
                     tol = 1e-3)$minimum
          }
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
cat("SIMULATION: ric2\n")



M           <- num_cores
p_feat      <- 500
n_train     <- 200
N_post      <- 3000
N_burn      <- 2000

sigma       <- 1
sim_par <- list(
  qs     = 5,
  qw     = 1, # starting value, will be changed by the functions
  cWe    = 5.0,
  cWf    = 1.0,
  cSe    = 0.5,
  cSf    = 2,
  tau0   = 1,
  hetero = FALSE
)


param_grid_for_CV <- expand.grid(
  # k = c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5)
  k = c(0.2, 0.25, 0.3)
)

qw_vals <- c(1, 4, 8, 12, 16, 20)
# registerDoParallel(cores = 6)
# cens_scales <- cens_scales_nonlinear <- compute_cens_scales_parallel(qw_vals, n_train, p_feat, sim_par, sigma, 0.35, 500)
cens_scales <- c(0.02357744, 0.02387869, 0.02674325, 0.02329564, 0.02420517, 0.02471282)

results <- run_simulations_over_qw(
  qw_values = qw_vals,
  cens_scales = cens_scales,
  M = M,
  n_train = n_train,
  p_feat = p_feat,
  sigma = sigma,
  par = sim_par,
  k = 0.1, 
  param_grid = param_grid_for_CV,
  N_post = N_post,
  N_burn = N_burn
)

info <- list(
  design = list(
    p_feat  = p_feat,
    n_train = n_train,
    sigma   = sigma,
    target_censoring = 0.35,
    qs      = sim_par$qs,
    qw_vals = qw_vals,
    cWe     = sim_par$cWe,
    cWf     = sim_par$cWf,
    cSe     = sim_par$cSe,
    cSf     = sim_par$cSf,
    tau0    = sim_par$tau0,
    hetero  = sim_par$hetero
  ),
  
  models = list(
    N_post = N_post,
    N_burn = N_burn,
    CHF_k_grid  = param_grid_for_CV$k
  ),
  
  computation = list(
    M = M,
    cores_used = num_cores
  )
)

combined_results <- list(results = results,
                         info = info)


# Define output file path 
# (! NAME MUST BE FILENAME_output.rds !)
output_file <- file.path(Sys.getenv("TMPDIR"), "ric2_output.rds")

# Print message
cat("Saving all settings results to:", output_file, "\n")

# Save the combined list
saveRDS(combined_results, file = output_file)

# Confirm successful save
cat("All results successfully saved in one file.\n")





























# 
# library(dplyr)
# library(tidyr)
# library(ggplot2)
# 
# res <- readRDS(
#   file.path(
#     "~/Library/CloudStorage/OneDrive-VrijeUniversiteitAmsterdam/Documents/GitHub/ShrinkageTrees/revision/RIC",
#     "ric2_output.rds"
#   )
# )
# 
# results <- res$results
# 
# summary_res <- results %>%
#   group_by(Method, qw) %>%
#   summarise(
#     RMSE     = mean(RMSE, na.rm = TRUE),
#     AbsBias  = mean(CATE_AbsBias, na.rm = TRUE),
#     Coverage = mean(CATE_coverage, na.rm = TRUE),
#     .groups = "drop"
#   )
# 
# summary_long <- summary_res %>%
#   pivot_longer(
#     cols = c(RMSE, AbsBias, Coverage),
#     names_to = "Metric",
#     values_to = "Value"
#   )
# 
# ggplot(summary_long, aes(x = qw, y = Value, color = Method)) +
#   geom_line(linewidth = 1) +
#   geom_point(size = 2) +
#   facet_wrap(~ Metric, scales = "free_y") +
#   labs(
#     x = "Number of weak confounders (qw)",
#     y = "ric2",
#     color = "Method"
#   ) +
#   theme_minimal() +
#   theme(
#     legend.position = "bottom",
#     strip.text = element_text(face = "bold")
#   )