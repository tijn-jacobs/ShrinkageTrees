library(ShrinkageTrees)
library(foreach)
library(doParallel)

source("evaluation_functions.R")


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
sigma <- sqrt(3)
s_prog <- 0.1
s_treat <- 0.05
N_post <- 3000
N_burn <- 2000

p_vals <- round(10^seq(log10(10), log10(10000), length.out = 25)) # Important to set number of p's!!!
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
  0.05735982, 0.05677074, 0.05680847, 0.05713565, 0.05596345,
  0.05650599, 0.05732633, 0.05640099, 0.05452684, 0.05314774,
  0.05360828, 0.05220717, 0.05048066, 0.05076701, 0.04674042,
  0.04446327, 0.04202325, 0.03979075, 0.03749388, 0.03554021,
  0.03449886, 0.03342433, 0.03182504, 0.03262220, 0.02984496
)

param_grid_for_CV <- expand.grid(
  k = c(0.1, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5)
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
output_file <- file.path(Sys.getenv("TMPDIR"), "quick0_output.rds")

# Print message
cat("Saving all settings results to:", output_file, "\n")

# Save the combined list
saveRDS(combined_results, file = output_file)

# Confirm successful save
cat("All results successfully saved in one file.\n")







medium <- readRDS("/Users/tijnjacobs/Library/CloudStorage/OneDrive-VrijeUniversiteitAmsterdam/Documents/GitHub/ShrinkageTrees/revision/quick0_output.rds")

medium <- medium$results_medium



library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

#==================================================
# Load medium results
#==================================================
medium <- medium %>%
  mutate(ATE_RMSE = ATE_bias^2)

metrics <- c(
  "sigma","ATE","ATE_bias","ATE_RMSE","ATE_coverage",
  "ATE_CI_Length","CATE_RPEHE","CATE_coverage",
  "CATE_CI_Length","RMSE","C_Index"
)

#==================================================
# Summaries
#==================================================
sum_medium <- medium %>%
  group_by(Method, p_feat) %>%
  summarise(across(all_of(metrics), mean, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(ATE_RMSE = sqrt(ATE_RMSE))

#==================================================
# Rename methods to uniform labels
#==================================================
rename_methods <- function(x) {
  recode(x,
         "IndivAFT"    = "AFT_BART",
         "BCF_AFT"     = "AFT_BCF",
         "DART_AFT"    = "AFT_DART",
         "S_BCF_AFT"   = "AFT_S_BCF",
         "CHF (k=0.1)" = "CHF",
         "CHF_CV"      = "CHF_CV")
}

sum_medium$Method <- rename_methods(sum_medium$Method)
medium$Method      <- rename_methods(medium$Method)

#==================================================
# Color palette (6 methods)
#==================================================
method_colors <- c(
  "AFT_BART"  = "#1f78b4",
  "AFT_BCF"   = "#33a02c",
  "AFT_DART"  = "#e31a1c",
  "AFT_S_BCF" = "#ff7f00",
  "CHF"       = "#6a3d9a",
  "CHF_CV"    = "#cab2d6"
)

linetype_values <- c(ATE = "dashed", CATE = "solid")
shape_values    <- c(ATE = 1, CATE = 0)

#==================================================
# Smoothing settings
#==================================================
smooth_fun <- geom_smooth(
  se = FALSE,
  method = "loess",
  span = 1.2,      # Adjust this to tune smoothness
  linewidth = 0.9
)

#==================================================
# 1) RMSE PLOT
#==================================================
long_rmse <- sum_medium %>%
  select(Method, p_feat, ATE_RMSE, CATE_RPEHE) %>%
  pivot_longer(c(ATE_RMSE, CATE_RPEHE),
               names_to = "Metric",
               values_to = "Value") %>%
  mutate(Metric = recode(Metric,
                         ATE_RMSE = "ATE",
                         CATE_RPEHE = "CATE"))

rmse_plot <- ggplot(long_rmse, aes(x = p_feat, y = Value,
                                   color = Method,
                                   linetype = Metric,
                                   shape = Metric)) +
  smooth_fun +
  scale_x_continuous(limits = c(0, 10000)) +
  scale_color_manual(values = method_colors) +
  scale_linetype_manual(values = linetype_values) +
  scale_shape_manual(values = shape_values) +
  labs(x = "Number of covariates", y = "RMSE", title = "RMSE") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    text = element_text(family = "Times New Roman"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

#==================================================
# 2) COVERAGE PLOT
#==================================================
sum_medium_cov <- medium %>%
  select(Method, p_feat, ATE_coverage, CATE_coverage) %>%
  group_by(Method, p_feat) %>%
  summarise(across(c(ATE_coverage, CATE_coverage), mean, na.rm = TRUE),
            .groups = "drop")

long_cov <- sum_medium_cov %>%
  pivot_longer(c(ATE_coverage, CATE_coverage),
               names_to = "Metric",
               values_to = "Value") %>%
  mutate(Metric = recode(Metric,
                         ATE_coverage = "ATE",
                         CATE_coverage = "CATE"))

coverage_plot <- ggplot(long_cov, aes(x = p_feat, y = Value,
                                      color = Method,
                                      linetype = Metric,
                                      shape = Metric)) +
  smooth_fun +
  scale_x_continuous(limits = c(0, 10000)) +
  scale_color_manual(values = method_colors) +
  scale_linetype_manual(values = linetype_values) +
  scale_shape_manual(values = shape_values) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Number of covariates", y = "Coverage", title = "Coverage") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    text = element_text(family = "Times New Roman"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

#==================================================
# 3) CI LENGTH PLOT
#==================================================
long_cilen <- sum_medium %>%
  select(Method, p_feat, ATE_CI_Length, CATE_CI_Length) %>%
  pivot_longer(c(ATE_CI_Length, CATE_CI_Length),
               names_to = "Metric",
               values_to = "Value") %>%
  mutate(Metric = recode(Metric,
                         ATE_CI_Length = "ATE",
                         CATE_CI_Length = "CATE"))

cilen_plot <- ggplot(long_cilen, aes(x = p_feat, y = Value,
                                     color = Method,
                                     linetype = Metric,
                                     shape = Metric)) +
  smooth_fun +
  scale_x_continuous(limits = c(0, 10000)) +
  scale_color_manual(values = method_colors) +
  scale_linetype_manual(values = linetype_values) +
  scale_shape_manual(values = shape_values) +
  labs(x = "Number of covariates", y = "CI length",
       title = "Credible Interval Length") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    text = element_text(family = "Times New Roman"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

#==================================================
# Combined plot (3 panels)
#==================================================
combined_plot <- rmse_plot + coverage_plot + cilen_plot + plot_layout(ncol = 3)

combined_plot



















library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

# ==============================
# X-axis scale choice (EDIT HERE)
# ==============================
x_scale_trans <- "identity"   # "log10", "log2", or "log"

#==================================================
# Prepare data
#==================================================
medium <- medium %>%
  mutate(ATE_RMSE = ATE_bias^2)

metrics <- c(
  "sigma","ATE","ATE_bias","ATE_RMSE","ATE_coverage",
  "ATE_CI_Length","CATE_RPEHE","CATE_coverage",
  "CATE_CI_Length","RMSE","C_Index"
)

#==================================================
# Summaries
#==================================================
sum_medium <- medium %>%
  group_by(Method, p_feat) %>%
  summarise(across(all_of(metrics), mean, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(ATE_RMSE = sqrt(ATE_RMSE))

#==================================================
# Rename methods
#==================================================
rename_methods <- function(x) {
  recode(x,
         "IndivAFT"    = "AFT_BART",
         "BCF_AFT"     = "AFT_BCF",
         "DART_AFT"    = "AFT_DART",
         "S_BCF_AFT"   = "AFT_S_BCF",
         "CHF (k=0.1)" = "CHF",
         "CHF_CV"      = "CHF_CV")
}

sum_medium$Method <- rename_methods(sum_medium$Method)
medium$Method     <- rename_methods(medium$Method)

#==================================================
# Color & style settings
#==================================================
method_colors <- c(
  "AFT_BART"  = "#1f78b4",
  "AFT_BCF"   = "#33a02c",
  "AFT_DART"  = "#e31a1c",
  "AFT_S_BCF" = "#ff7f00",
  "CHF"       = "#6a3d9a",
  "CHF_CV"    = "#cab2d6"
)

linetype_values <- c(ATE = "dashed", CATE = "solid")
shape_values    <- c(ATE = 1, CATE = 0)

#==================================================
# Smoothing
#==================================================
smooth_fun <- geom_smooth(
  se = FALSE,
  method = "loess",
  span = 1.2,
  linewidth = 0.9
)

#==================================================
# 1) RMSE plot
#==================================================
long_rmse <- sum_medium %>%
  select(Method, p_feat, ATE_RMSE, CATE_RPEHE) %>%
  pivot_longer(c(ATE_RMSE, CATE_RPEHE),
               names_to = "Metric",
               values_to = "Value") %>%
  mutate(Metric = recode(Metric,
                         ATE_RMSE = "ATE",
                         CATE_RPEHE = "CATE"))

rmse_plot <- ggplot(long_rmse,
                    aes(x = p_feat, y = Value,
                        color = Method,
                        linetype = Metric,
                        shape = Metric)) +
  smooth_fun +
  scale_x_continuous(
    trans  = x_scale_trans,
    breaks = sort(unique(long_rmse$p_feat))
  ) +
  scale_color_manual(values = method_colors) +
  scale_linetype_manual(values = linetype_values) +
  scale_shape_manual(values = shape_values) +
  labs(x = "Number of covariates", y = "RMSE", title = "RMSE") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    text = element_text(family = "Times New Roman"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

#==================================================
# 2) Coverage plot
#==================================================
sum_medium_cov <- medium %>%
  select(Method, p_feat, ATE_coverage, CATE_coverage) %>%
  group_by(Method, p_feat) %>%
  summarise(across(c(ATE_coverage, CATE_coverage), mean, na.rm = TRUE),
            .groups = "drop")

long_cov <- sum_medium_cov %>%
  pivot_longer(c(ATE_coverage, CATE_coverage),
               names_to = "Metric",
               values_to = "Value") %>%
  mutate(Metric = recode(Metric,
                         ATE_coverage = "ATE",
                         CATE_coverage = "CATE"))

coverage_plot <- ggplot(long_cov,
                        aes(x = p_feat, y = Value,
                            color = Method,
                            linetype = Metric,
                            shape = Metric)) +
  smooth_fun +
  scale_x_continuous(
    trans  = x_scale_trans,
    breaks = sort(unique(long_cov$p_feat))
  ) +
  scale_color_manual(values = method_colors) +
  scale_linetype_manual(values = linetype_values) +
  scale_shape_manual(values = shape_values) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Number of covariates", y = "Coverage", title = "Coverage") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    text = element_text(family = "Times New Roman"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

#==================================================
# 3) CI length plot
#==================================================
long_cilen <- sum_medium %>%
  select(Method, p_feat, ATE_CI_Length, CATE_CI_Length) %>%
  pivot_longer(c(ATE_CI_Length, CATE_CI_Length),
               names_to = "Metric",
               values_to = "Value") %>%
  mutate(Metric = recode(Metric,
                         ATE_CI_Length = "ATE",
                         CATE_CI_Length = "CATE"))

cilen_plot <- ggplot(long_cilen,
                     aes(x = p_feat, y = Value,
                         color = Method,
                         linetype = Metric,
                         shape = Metric)) +
  smooth_fun +
  scale_x_continuous(
    trans  = x_scale_trans,
    breaks = sort(unique(long_cilen$p_feat))
  ) +
  scale_color_manual(values = method_colors) +
  scale_linetype_manual(values = linetype_values) +
  scale_shape_manual(values = shape_values) +
  labs(x = "Number of covariates", y = "CI length",
       title = "Credible Interval Length") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    text = element_text(family = "Times New Roman"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

#==================================================
# Combined plot
#==================================================
combined_plot <- rmse_plot + coverage_plot + cilen_plot +
  plot_layout(ncol = 3)

combined_plot
