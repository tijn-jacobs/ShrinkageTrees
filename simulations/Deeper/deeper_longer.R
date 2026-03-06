### Simulation deeper Friedman ###

library(ShrinkageTrees)
library(foreach)
library(doParallel)

source("evaluation_functions.R")


data_gen <- function(n_train, p_feat, sigma, s_prog, s_treat = NULL, cens_scale = 1, linear = NULL) {
  X_train <- matrix(runif(n_train * p_feat), n_train, p_feat)
  
  propensity <- pnorm(-0.5 + 0.4*X_train[,1])
  treatment <- rbinom(n_train, 1, propensity)
  
  beta_prog <- rnorm(p_feat, 0, 1) * rbinom(p_feat, 1, s_prog)
  mu <- beta_prog %*% t(X_train)
  
  tau <- 10 * sin(pi * X_train[, 1] * X_train[, 2]) +
    20 * (X_train[, 3] - 0.5)^2 +
    10 * X_train[, 4] +
    5  * X_train[, 5]
  
  true_ate <- 10 * 0.5246745 + 20 * 1/12 + 10 * 1/2 + 5 * 1/2
  
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
                            s_prog = 0.1,
                            target_cens,
                            M = 1000,
                            verbose = FALSE) {

  # Inner function: return absolute difference between actual and target censoring
  estimate_cens_diff <- function(cens_scale) {
    rates <- replicate(M, {
      dt <- data_gen(n_train, p, sigma, s_prog, s_treat = NULL, cens_scale, linear = NULL)
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
sigma <- sqrt(3)
s_prog <- 0.1
N_post <- 3
N_burn <- 2

p_vals <- round(10^seq(log10(10), log10(2500), length.out = 20)) # Important to set number of p's!!!
# cens_scales_medium <- sapply(p_vals, function(p) find_cens_scale(p, n = 200, target_cens = 0.35, M = 500))
cens_scales_medium <- c(0.06226613, 0.06984927, 0.06043309, 0.07018260, 0.06084431, 0.06043760, 0.05894593, 0.05958603, 0.05929964, 0.05837365, 0.05754537, 0.05566082, 0.05437630, 0.05348404, 0.05110670, 0.04875770, 0.04600227, 0.04402867, 0.04074547, 0.03743365)

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
  param_grid = param_grid_for_CV,
  k = 0.1, # Redundant
  N_post = N_post,
  N_burn = N_burn
)

combined_results <- list(results_medium = results_medium)


# Define output file path 
# (! NAME MUST BE FILENAME_output.rds !)
output_file <- file.path(Sys.getenv("TMPDIR"), "deeper_longer_output.rds")

# Print message
cat("Saving all settings results to:", output_file, "\n")

# Save the combined list
saveRDS(combined_results, file = output_file)

# Confirm successful save
cat("All results successfully saved in one file.\n")





# 
# 
# medium_longer <- readRDS("/Users/tijnjacobs/Library/CloudStorage/OneDrive-VrijeUniversiteitAmsterdam/Documents/GitHub/ShrinkageTrees/simulations revision/deeper_longer_output.rds")
# medium_all <- readRDS("/Users/tijnjacobs/Library/CloudStorage/OneDrive-VrijeUniversiteitAmsterdam/Documents/GitHub/ShrinkageTrees/simulations revision/deeper_longer_output.rds")
# medium <- medium$results_medium
# 
# 
# 
# 
# library(dplyr)
# library(tidyr)
# library(ggplot2)
# library(patchwork)
# 
# ### --- Prepare medium dataset --- ###
# 
# medium <- medium %>%
#   mutate(ATE_RMSE = ATE_bias^2)
# 
# metrics <- c(
#   "sigma","ATE","ATE_bias","ATE_RMSE","ATE_coverage",
#   "ATE_CI_Length","CATE_RPEHE","CATE_coverage",
#   "CATE_CI_Length","RMSE","C_Index"
# )
# 
# # Summaries
# sum_medium <- medium %>%
#   group_by(Method, p_feat) %>%
#   summarise(across(all_of(metrics), ~ mean(.x, na.rm = TRUE)),
#             .groups = "drop") %>%
#   mutate(ATE_RMSE = sqrt(ATE_RMSE))
# 
# ### --- Rename methods uniformly --- ###
# 
# rename_methods <- function(x) {
#   recode(x,
#          "IndivAFT"    = "AFT_BART",
#          "BCF_AFT"     = "AFT_BCF",
#          "DART_AFT"    = "AFT_DART",
#          "S_BCF_AFT"   = "AFT_S_BCF",
#          "CHF (k=0.1)" = "CHF",
#          "CHF_CV"      = "CHF_CV")
# }
# 
# sum_medium$Method <- rename_methods(sum_medium$Method)
# 
# ### --- Colors --- ###
# 
# method_colors <- c(
#   "AFT_BART"  = "#1f78b4",
#   "AFT_BCF"   = "#33a02c",
#   "AFT_DART"  = "#e31a1c",
#   "AFT_S_BCF" = "#ff7f00",
#   "CHF"       = "#6a3d9a",
#   "CHF_CV"    = "#cab2d6"
# )
# 
# linetype_values <- c(ATE = "dashed", CATE = "solid")
# shape_values    <- c(ATE = 1, CATE = 0)
# 
# # ===============================================================
# #                     1) RMSE PLOT
# # ===============================================================
# 
# long_rmse <- sum_medium %>%
#   dplyr::select(Method, p_feat, ATE_RMSE, CATE_RPEHE) %>%
#   tidyr::pivot_longer(
#     c(ATE_RMSE, CATE_RPEHE),
#     names_to = "Metric",
#     values_to = "Value"
#   ) %>%
#   mutate(Metric = recode(Metric,
#                          ATE_RMSE = "ATE",
#                          CATE_RPEHE = "CATE"))
# 
# rmse_plot <- ggplot(long_rmse,
#                     aes(x = p_feat, y = Value,
#                         color = Method, linetype = Metric,
#                         shape = Metric)) +
#   geom_smooth(se = FALSE, linewidth = 0.8) +
#   scale_x_continuous(limits = c(0, 2500)) +
#   scale_color_manual(values = method_colors) +
#   scale_linetype_manual(values = linetype_values) +
#   scale_shape_manual(values = shape_values) +
#   labs(
#     x = "Number of covariates",
#     y = "RMSE",
#     title = "RMSE"
#   ) +
#   theme_minimal(base_size = 14) +
#   theme(
#     legend.position = "none",
#     text = element_text(family = "Times New Roman"),
#     axis.text.x = element_text(angle = 45, hjust = 1)
#   )
# 
# # ===============================================================
# #                     2) COVERAGE PLOT
# # ===============================================================
# 
# sum_medium_cov <- medium %>%
#   dplyr::select(Method, p_feat, ATE_coverage, CATE_coverage) %>%
#   group_by(Method, p_feat) %>%
#   summarise(across(c(ATE_coverage, CATE_coverage), mean, na.rm = TRUE),
#             .groups = "drop")
# 
# sum_medium_cov$Method <- rename_methods(sum_medium_cov$Method)
# 
# long_cov <- sum_medium_cov %>%
#   tidyr::pivot_longer(
#     c(ATE_coverage, CATE_coverage),
#     names_to = "Metric",
#     values_to = "Value"
#   ) %>%
#   mutate(Metric = recode(Metric,
#                          ATE_coverage = "ATE",
#                          CATE_coverage = "CATE"))
# 
# coverage_plot <- ggplot(long_cov,
#                         aes(x = p_feat, y = Value,
#                             color = Method, linetype = Metric,
#                             shape = Metric)) +
#   geom_smooth(se = FALSE, linewidth = 0.8) +
#   scale_x_continuous(limits = c(0, 2500)) +
#   scale_color_manual(values = method_colors) +
#   scale_linetype_manual(values = linetype_values) +
#   scale_shape_manual(values = shape_values) +
#   coord_cartesian(ylim = c(0, 1)) +
#   labs(
#     x = "Number of covariates",
#     y = "Coverage",
#     title = "Coverage"
#   ) +
#   theme_minimal(base_size = 14) +
#   theme(
#     legend.position = "none",
#     text = element_text(family = "Times New Roman"),
#     axis.text.x = element_text(angle = 45, hjust = 1)
#   )
# 
# # ===============================================================
# #                     3) CI LENGTH PLOT
# # ===============================================================
# 
# long_cilen <- sum_medium %>%
#   dplyr::select(Method, p_feat, ATE_CI_Length, CATE_CI_Length) %>%
#   tidyr::pivot_longer(
#     c(ATE_CI_Length, CATE_CI_Length),
#     names_to = "Metric",
#     values_to = "Value"
#   ) %>%
#   mutate(Metric = recode(Metric,
#                          ATE_CI_Length = "ATE",
#                          CATE_CI_Length = "CATE"))
# 
# cilen_plot <- ggplot(long_cilen,
#                      aes(x = p_feat, y = Value,
#                          color = Method, linetype = Metric,
#                          shape = Metric)) +
#   geom_smooth(se = FALSE, linewidth = 0.8) +
#   scale_x_continuous(limits = c(0, 2500)) +
#   scale_color_manual(values = method_colors) +
#   scale_linetype_manual(values = linetype_values) +
#   scale_shape_manual(values = shape_values) +
#   labs(
#     x = "Number of covariates",
#     y = "CI length",
#     title = "Credible Interval Length"
#   ) +
#   theme_minimal(base_size = 14) +
#   theme(
#     legend.position = "right",
#     text = element_text(family = "Times New Roman"),
#     axis.text.x = element_text(angle = 45, hjust = 1)
#   )
# 
# # ===============================================================
# #                     COMBINE ALL 3 PLOTS
# # ===============================================================
# 
# combined_plot <- rmse_plot + coverage_plot + cilen_plot + plot_layout(ncol = 3)
# 
# combined_plot
# 
