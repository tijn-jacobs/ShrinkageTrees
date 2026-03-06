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

M <- 1000
n_train <- 100
sigma <- sqrt(3)
s_prog <- 0.1
s_treat <- 0.05
N_post <- 3000
N_burn <- 2000

make_p_vals <- function(p_min = 10, p_max = 5000,
                        n_log = 12, n_high_linear = 18,
                        high_from = 1000) {
  log_part  <- round(10^seq(log10(p_min), log10(p_max), length.out = n_log))
  high_part <- round(seq(high_from, p_max, length.out = n_high_linear))
  sort(unique(pmin(p_max, pmax(p_min, c(log_part, high_part)))))
}

p_vals <- make_p_vals()

# cens_scales_medium <- foreach(p = p_vals, .combine = c) %dopar% {
#   find_cens_scale(
#     p            = p,
#     n_train      = n_train,
#     sigma        = sigma,
#     s_prog       = s_prog,
#     s_treat      = s_treat,
#     target_cens  = 0.35,
#     M            = 500,
#     verbose      = FALSE
#   )
# }

cens_scales_medium <- c(
  0.06066867, 0.06108191, 0.06201221, 0.05897433, 0.05877816,
  0.05730583, 0.05525392, 0.05008671, 0.04652164, 0.04466800,
  0.04355119, 0.04033435, 0.04478236, 0.03916834, 0.03886713,
  0.03710918, 0.03673590, 0.03590992, 0.03396221, 0.03400513,
  0.03424924, 0.03320158, 0.02696475, 0.03384381, 0.03176677,
  0.03094370, 0.03139459, 0.03060235, 0.03039377
)


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
  k = 0.1, # Redundant
  N_post = N_post,
  N_burn = N_burn
)

combined_results <- list(results_medium = results_medium)


# Define output file path 
# (! NAME MUST BE FILENAME_output.rds !)
output_file <- file.path(Sys.getenv("TMPDIR"), "revision_deeper_output.rds")

# Print message
cat("Saving all settings results to:", output_file, "\n")

# Save the combined list
saveRDS(combined_results, file = output_file)

# Confirm successful save
cat("All results successfully saved in one file.\n")






# 
# medium <- readRDS("/Users/tijnjacobs/Library/CloudStorage/OneDrive-VrijeUniversiteitAmsterdam/Documents/GitHub/ShrinkageTrees/revision/deeper/revision_deeper_output.rds")
# 
# medium <- medium$results_medium
# 
# library(dplyr)
# library(tidyr)
# library(ggplot2)
# library(patchwork)
# 
# #==================================================
# # METHOD SELECTION (EDIT HERE)
# #==================================================
# # methods_to_plot <- c(
# #   "CHF_CV", "CHF0.1", "CHF0.2", "CHF0.25", "CHF0.3", "CHF0.5",
# #   "AFT_BART", "AFT_DART", "AFT_BCF", "AFT_BCF_shrink"
# # )
# 
# methods_to_plot <- c(
#   "CHF_CV", "CHF0.1", "CHF0.2", "CHF0.25", "CHF0.3", "CHF0.5",
#   "AFT_BCF", "AFT_BCF_shrink"
# )
# 
# #==================================================
# # Filter data
# #==================================================
# medium_plot <- medium %>%
#   filter(Method %in% methods_to_plot)
# 
# #==================================================
# # Summaries (CATE metrics only)
# #==================================================
# sum_medium <- medium_plot %>%
#   group_by(Method, p_feat) %>%
#   summarise(
#     CATE_RPEHE    = mean(CATE_RPEHE, na.rm = TRUE),
#     CATE_coverage = mean(CATE_coverage, na.rm = TRUE),
#     CATE_CI_Length = mean(CATE_CI_Length, na.rm = TRUE),
#     .groups = "drop"
#   )
# 
# #==================================================
# # Color palette (one color per method)
# #==================================================
# method_colors <- c(
#   "CHF_CV"         = "#000000",  # black (anchor / reference)
#   "CHF0.1"         = "#4E79A7",  # blue
#   "CHF0.2"         = "#59A14F",  # green
#   "CHF0.25"        = "#F28E2B",  # orange
#   "CHF0.3"         = "#E15759",  # red
#   "CHF0.5"         = "#B07AA1",  # purple
#   
#   "AFT_BART"       = "#76B7B2",  # teal
#   "AFT_DART"       = "#EDC948",  # yellow
#   "AFT_BCF"        = "#9C755F",  # brown
#   "AFT_BCF_shrink" = "#FF9DA7"   # pink
# )
# 
# # Keep only needed colors (important when subsetting)
# method_colors <- method_colors[names(method_colors) %in% methods_to_plot]
# 
# #==================================================
# # Smoothing
# #==================================================
# smooth_fun <- geom_smooth(
#   se = FALSE,
#   method = "loess",
#   span = 1.1,
#   linewidth = 0.9
# )
# 
# smooth_fun <- geom_smooth(
#   se = FALSE,
#   method = "gam",
#   formula = y ~ s(x, bs = "cs"),
#   linewidth = 0.9
# )
# 
# #==================================================
# # 1) CATE RMSE (RPEHE)
# #==================================================
# rmse_plot <- ggplot(
#   sum_medium,
#   aes(x = p_feat, y = CATE_RPEHE, color = Method)
# ) +
#   smooth_fun +
#   scale_x_continuous(
#     trans  = "log10",
#     breaks = c(50, 500, 5000),
#     limits = c(10, 5000)
#   ) +
#   scale_color_manual(values = method_colors) +
#   labs(x = "Number of covariates", y = "CATE RMSE", title = "CATE RMSE") +
#   theme_minimal(base_size = 14) +
#   theme(
#     legend.position = "none",
#     text = element_text(family = "Times New Roman"),
#     axis.text.x = element_text(angle = 45, hjust = 1)
#   )
# 
# #==================================================
# # 2) CATE Coverage
# #==================================================
# coverage_plot <- ggplot(
#   sum_medium,
#   aes(x = p_feat, y = CATE_coverage, color = Method)
# ) +
#   smooth_fun +
#   geom_hline(yintercept = 0.95, linetype = "dashed", color = "black") +
#   scale_x_continuous(
#     trans  = "log10",
#     breaks = c(50, 500, 5000),
#     limits = c(10, 5000)
#   ) +
#   scale_color_manual(values = method_colors) +
#   coord_cartesian(ylim = c(0, 1)) +
#   labs(x = "Number of covariates", y = "Coverage", title = "CATE Coverage") +
#   theme_minimal(base_size = 14) +
#   theme(
#     legend.position = "none",
#     text = element_text(family = "Times New Roman"),
#     axis.text.x = element_text(angle = 45, hjust = 1)
#   )
# 
# #==================================================
# # 3) CATE CI Length
# #==================================================
# cilen_plot <- ggplot(
#   sum_medium,
#   aes(x = p_feat, y = CATE_CI_Length, color = Method)
# ) +
#   smooth_fun +
#   scale_x_continuous(
#     trans  = "log10",
#     breaks = c(50, 500, 5000),
#     limits = c(10, 5000)
#   ) +
#   scale_color_manual(values = method_colors) +
#   labs(x = "Number of covariates", y = "CI length",
#        title = "CATE Credible Interval Length") +
#   theme_minimal(base_size = 14) +
#   theme(
#     legend.position = "right",
#     text = element_text(family = "Times New Roman"),
#     axis.text.x = element_text(angle = 45, hjust = 1)
#   )
# 
# #==================================================
# # Combined plot
# #==================================================
# combined_plot <- rmse_plot + coverage_plot + cilen_plot +
#   plot_layout(ncol = 3)
# 
# combined_plot
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# library(dplyr)
# library(tidyr)
# library(ggplot2)
# library(patchwork)
# 
# # ==============================
# # X-axis scale choice (EDIT HERE)
# # ==============================
# x_scale_trans <- "identity"   # "log10", "log2", or "log"
# 
# #==================================================
# # Prepare data
# #==================================================
# medium <- medium %>%
#   mutate(ATE_RMSE = ATE_bias^2)
# 
# metrics <- c(
#   "sigma","ATE","ATE_bias","ATE_RMSE","ATE_coverage",
#   "ATE_CI_Length","CATE_RPEHE","CATE_coverage",
#   "CATE_CI_Length","RMSE","C_Index"
# )
# 
# #==================================================
# # Summaries
# #==================================================
# sum_medium <- medium %>%
#   group_by(Method, p_feat) %>%
#   summarise(across(all_of(metrics), mean, na.rm = TRUE),
#             .groups = "drop") %>%
#   mutate(ATE_RMSE = sqrt(ATE_RMSE))
# 
# #==================================================
# # Rename methods
# #==================================================
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
# medium$Method     <- rename_methods(medium$Method)
# 
# #==================================================
# # Color & style settings
# #==================================================
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
# #==================================================
# # Smoothing
# #==================================================
# smooth_fun <- geom_smooth(
#   se = FALSE,
#   method = "loess",
#   span = 1.2,
#   linewidth = 0.9
# )
# 
# #==================================================
# # 1) RMSE plot
# #==================================================
# long_rmse <- sum_medium %>%
#   select(Method, p_feat, ATE_RMSE, CATE_RPEHE) %>%
#   pivot_longer(c(ATE_RMSE, CATE_RPEHE),
#                names_to = "Metric",
#                values_to = "Value") %>%
#   mutate(Metric = recode(Metric,
#                          ATE_RMSE = "ATE",
#                          CATE_RPEHE = "CATE"))
# 
# rmse_plot <- ggplot(long_rmse,
#                     aes(x = p_feat, y = Value,
#                         color = Method,
#                         linetype = Metric,
#                         shape = Metric)) +
#   smooth_fun +
#   scale_x_continuous(
#     trans  = x_scale_trans,
#     breaks = sort(unique(long_rmse$p_feat))
#   ) +
#   scale_color_manual(values = method_colors) +
#   scale_linetype_manual(values = linetype_values) +
#   scale_shape_manual(values = shape_values) +
#   labs(x = "Number of covariates", y = "RMSE", title = "RMSE") +
#   theme_minimal(base_size = 14) +
#   theme(
#     legend.position = "none",
#     text = element_text(family = "Times New Roman"),
#     axis.text.x = element_text(angle = 45, hjust = 1)
#   )
# 
# #==================================================
# # 2) Coverage plot
# #==================================================
# sum_medium_cov <- medium %>%
#   select(Method, p_feat, ATE_coverage, CATE_coverage) %>%
#   group_by(Method, p_feat) %>%
#   summarise(across(c(ATE_coverage, CATE_coverage), mean, na.rm = TRUE),
#             .groups = "drop")
# 
# long_cov <- sum_medium_cov %>%
#   pivot_longer(c(ATE_coverage, CATE_coverage),
#                names_to = "Metric",
#                values_to = "Value") %>%
#   mutate(Metric = recode(Metric,
#                          ATE_coverage = "ATE",
#                          CATE_coverage = "CATE"))
# 
# coverage_plot <- ggplot(long_cov,
#                         aes(x = p_feat, y = Value,
#                             color = Method,
#                             linetype = Metric,
#                             shape = Metric)) +
#   smooth_fun +
#   scale_x_continuous(
#     trans  = x_scale_trans,
#     breaks = sort(unique(long_cov$p_feat))
#   ) +
#   scale_color_manual(values = method_colors) +
#   scale_linetype_manual(values = linetype_values) +
#   scale_shape_manual(values = shape_values) +
#   coord_cartesian(ylim = c(0, 1)) +
#   labs(x = "Number of covariates", y = "Coverage", title = "Coverage") +
#   theme_minimal(base_size = 14) +
#   theme(
#     legend.position = "none",
#     text = element_text(family = "Times New Roman"),
#     axis.text.x = element_text(angle = 45, hjust = 1)
#   )
# 
# #==================================================
# # 3) CI length plot
# #==================================================
# long_cilen <- sum_medium %>%
#   select(Method, p_feat, ATE_CI_Length, CATE_CI_Length) %>%
#   pivot_longer(c(ATE_CI_Length, CATE_CI_Length),
#                names_to = "Metric",
#                values_to = "Value") %>%
#   mutate(Metric = recode(Metric,
#                          ATE_CI_Length = "ATE",
#                          CATE_CI_Length = "CATE"))
# 
# cilen_plot <- ggplot(long_cilen,
#                      aes(x = p_feat, y = Value,
#                          color = Method,
#                          linetype = Metric,
#                          shape = Metric)) +
#   smooth_fun +
#   scale_x_continuous(
#     trans  = x_scale_trans,
#     breaks = sort(unique(long_cilen$p_feat))
#   ) +
#   scale_color_manual(values = method_colors) +
#   scale_linetype_manual(values = linetype_values) +
#   scale_shape_manual(values = shape_values) +
#   labs(x = "Number of covariates", y = "CI length",
#        title = "Credible Interval Length") +
#   theme_minimal(base_size = 14) +
#   theme(
#     legend.position = "right",
#     text = element_text(family = "Times New Roman"),
#     axis.text.x = element_text(angle = 45, hjust = 1)
#   )
# 
# #==================================================
# # Combined plot
# #==================================================
# combined_plot <- rmse_plot + coverage_plot + cilen_plot +
#   plot_layout(ncol = 3)
# 
# combined_plot
