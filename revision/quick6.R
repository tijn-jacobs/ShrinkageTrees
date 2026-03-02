library(ShrinkageTrees)
library(foreach)
library(doParallel)

source("evaluation_functions.R")


data_gen <- function(n_train, p_feat, sigma,
                     s_prog, s_treat, cens_scale,
                     linear = NULL) {
  
  X_train <- matrix(runif(n_train * p_feat), n_train, p_feat)
  
  # Propensity & treatment (unchanged)
  propensity <- pnorm(-0.5 + 0.4 * X_train[, 1])
  treatment  <- rbinom(n_train, 1, propensity)
  
  # Prognostic component (unchanged)
  beta_prog <- rnorm(p_feat, 0, 1) * sort(rbinom(p_feat, 1, s_prog), TRUE)
  mu <- as.numeric(beta_prog %*% t(X_train))
  
  # --- Sparse threshold-based CATE ---
  m <- max(1, round(p_feat * s_treat))   # number of active covariates
  idx <- sample.int(p_feat, m)
  
  beta <- rnorm(m, mean = 0, sd = 2)     # moderate effect sizes
  thr  <- runif(m, 0.6, 0.9)              # thresholds make effects rare
  
  tau <- rep(0, n_train)
  true_ate_raw <- 0
  
  for (j in seq_len(m)) {
    active <- X_train[, idx[j]] > thr[j]
    tau[active] <- tau[active] + beta[j]
    
    # Analytic contribution: P(X > t) = 1 - t
    true_ate_raw <- true_ate_raw + beta[j] * (1 - thr[j])
  }
  
  # Event times + noise
  true_event_times <- mu + (treatment - 0.5) * tau
  uncensored_event_times <- true_event_times + rnorm(n_train, 0, sigma)
  
  sd_un <- sd(uncensored_event_times)
  uncensored_event_times <- uncensored_event_times / sd_un
  true_event_times       <- true_event_times / sd_un
  
  # Exponential censoring on log scale (unchanged)
  C <- log(rexp(n_train, cens_scale)) + min(uncensored_event_times)
  follow_up <- pmin(uncensored_event_times, C)
  status    <- as.numeric(uncensored_event_times <= C)
  
  list(
    X_train                = X_train,
    treatment              = treatment,
    propensity             = propensity,
    follow_up              = as.numeric(follow_up),
    status                 = status,
    true_event_times       = as.numeric(true_event_times),
    uncensored_event_times = as.numeric(uncensored_event_times),
    true_cate              = as.vector(tau) / sd_un,
    true_ate               = true_ate_raw / sd_un,
    sample_ate             = mean(tau) / sd_un,
    obs_sigma              = sigma / sd_un
  )
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
cat("SIMULATION: quick6\n")

M <- num_cores
n_train <- 100
sigma <- 1
s_prog <- 0.1
s_treat <- 0.1
N_post <- 3000
N_burn <- 2000

p_vals <- round(seq(10, 5000, length.out = 25)) # Important to set number of p's!!!
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
  0.02662624, 0.02854889, 0.02790231, 0.02763141, 0.02729808,
  0.02815583, 0.02867255, 0.02729808, 0.02729808, 0.02786768,
  0.02696475, 0.02815758, 0.02729808, 0.02899702, 0.02798749,
  0.02798394, 0.02729808, 0.02903784, 0.02910870, 0.02971435,
  0.02878078, 0.02862856, 0.02729808, 0.02820005, 0.02886250
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
output_file <- file.path(Sys.getenv("TMPDIR"), "quick6_output.rds")

# Print message
cat("Saving all settings results to:", output_file, "\n")

# Save the combined list
saveRDS(combined_results, file = output_file)

# Confirm successful save
cat("All results successfully saved in one file.\n")
















medium <- readRDS("/Users/tijnjacobs/Library/CloudStorage/OneDrive-VrijeUniversiteitAmsterdam/Documents/GitHub/ShrinkageTrees/simulations revision/quick6_output.rds")
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
  scale_x_continuous(limits = c(0, 5000)) +
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
  scale_x_continuous(limits = c(0, 5000)) +
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
  scale_x_continuous(limits = c(0, 5000)) +
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




