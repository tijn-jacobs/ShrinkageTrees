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
cat("SIMULATION: ric5\n")


M           <- 1000
p_feat      <- 500
n_train     <- 100
N_post      <- 5000
N_burn      <- 2500

sigma       <- 1
sim_par <- list(
  qs     = 5,
  qw     = 1, # starting value, will be changed by the functions
  cWe    = 2.0,
  cWf    = 1.0,
  cSe    = 2,
  cSf    = 4,
  tau0   = 1,
  hetero = TRUE
)


param_grid_for_CV <- expand.grid(
  # k = c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5)
  k = c(0.2, 0.25, 0.3)
)

qw_vals <- c(1:20)
# registerDoParallel(cores = 5)
# cens_scales <- cens_scales_nonlinear <- compute_cens_scales_parallel(qw_vals, n_train, p_feat, sim_par, sigma, 0.35, 200)
cens_scales <- c(
  0.03401976, 0.03265074, 0.03212704, 0.03418468, 0.03148633,
  0.03193255, 0.02607615, 0.03506591, 0.03304686, 0.03266444,
  0.03398019, 0.03316470, 0.03198963, 0.03377872, 0.03193707,
  0.03194194, 0.03118730, 0.03164874, 0.03339481, 0.03204421
)

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
output_file <- file.path(Sys.getenv("TMPDIR"), "ric5_output.rds")

# Print message
cat("Saving all settings results to:", output_file, "\n")

# Save the combined list
saveRDS(combined_results, file = output_file)

# Confirm successful save
cat("All results successfully saved in one file.\n")

























# 
# 
# 
# #==================================================
# # Load Libraries
# #==================================================
# library(dplyr)
# library(tidyr)
# library(ggplot2)
# library(ggh4x) # For facetted_pos_scales
# 
# #==================================================
# # Load Data
# #==================================================
# base_dir <- "~/Library/CloudStorage/OneDrive-VrijeUniversiteitAmsterdam/Documents/GitHub/ShrinkageTrees/revision/RIC"
# res <- readRDS(file.path(base_dir, "ric5_output.rds"))
# results <- res$results
# 
# 
# results_clean <- results %>%
#   filter(Method != "CHF0.15") %>% # Omit CHF0.2
#   mutate(Method = recode(Method,
#                          "CHF0.2"          = "Causal Horseshoe Forest",
#                          "AFT_BCF"          = "AFT-BCF",
#                          "AFT_ShrinkageBCF" = "AFT-Shrinkage BCF"
#   ))
# 
# # 2. Summarise
# summary_res <- results_clean %>%
#   group_by(Method, qw) %>%
#   summarise(
#     # We focus on RMSE and Coverage as requested
#     RMSE  = mean(CATE_RPEHE, na.rm = TRUE),
#     Coverage = mean(CATE_coverage, na.rm = TRUE),
#     .groups = "drop"
#   )
# 
# # 3. Pivot to Long Format
# summary_long <- summary_res %>%
#   pivot_longer(
#     cols = c(RMSE, Coverage),
#     names_to = "Metric",
#     values_to = "Value"
#   )
# 
# # 4. Force Panel Order (RMSE -> Coverage)
# summary_long$Metric <- factor(
#   summary_long$Metric,
#   levels = c("RMSE", "Coverage")
# )
# 
# #==================================================
# # Styling & Theme (Matching Reference)
# #==================================================
# # Define colors
# method_colors <- c(
#   "Causal Horseshoe Forest" = "#042E78",
#   "AFT-BCF"                = "#FFB74D",
#   "AFT-Shrinkage BCF"       = "#DD6E0F"
# )
# 
# # Define Theme
# theme_chf <- theme_bw(base_size = 16) +
#   theme(
#     text = element_text(family = "Times New Roman"),
#     axis.title.x = element_text(face = "italic"), # Italic x-axis title
#     axis.title.y = element_text(face = "plain"),
#     axis.text.x  = element_text(),
#     axis.text.y  = element_text(),
#     panel.border = element_rect(linewidth = 0.6),
#     legend.position = "right", # Bottom legend for this plot
#     legend.title = element_text(face = "plain", size = 16),
#     legend.key.width = unit(35, "pt"),
#     legend.text = element_text(size = 16),
#     strip.text = element_text(face = "bold", size = 20), # Bold panel titles
#     strip.background = element_blank()
#   )
# 
# #==================================================
# # Plotting
# #==================================================
# ric_plot <- ggplot(summary_long, aes(x = qw, y = Value, color = Method)) +
#   # Loess smoothing (matching reference style)
#   geom_smooth(
#     method = "loess",
#     se = FALSE,
#     span = 1.0,
#     linewidth = 1.2
#   ) +
#   # Points
#   # geom_point(size = 2.5, alpha = 0.9) +
# 
#   # Dashed line for Coverage only
#   geom_hline(
#     data = data.frame(Metric = factor("Coverage", levels = c("RMSE", "Coverage"))),
#     aes(yintercept = 0.95),
#     linetype = "dashed",
#     linewidth = 0.6,
#     color = "black"
#   ) +
# 
#   # Faceting
#   facet_wrap(~ Metric, scales = "free_y", nrow = 1) +
# 
#   scale_x_continuous(breaks = c(1, 10, 20)) + # Set X breaks here
#   facetted_pos_scales(
#     y = list(
#       Metric == "RMSE" ~ scale_y_continuous(
#         limits = c(0.45, 1.1),
#         breaks = c(0.6, 0.8, 1.0) # Manual Y breaks for RMSE
#       ),
#       Metric == "Coverage" ~ scale_y_continuous(
#         limits = c(0.0, 1.0),
#         breaks = c(0, 0.25, 0.5, 0.75, 1.0) # Manual Y breaks for Coverage
#       )
#     )
#   )+
# 
#   # Colors
#   scale_color_manual(values = method_colors) +
# 
#   # Labels with Math Expression for X-axis
#   labs(
#     x = "Number of weak confounders",
#     y = NULL,
#     color = "Method"
#   ) +
# 
#   # Apply Theme
#   theme_chf +
# 
#   # Thicker line in legend
#   guides(
#     color = guide_legend(
#       override.aes = list(linewidth = 2)
#     )
#   )
# 
# #==================================================
# # Save
# #==================================================
# ggsave(
#   filename = file.path(base_dir, "revision_ric.png"),
#   plot = ric_plot,
#   width = 1.1 * 250,
#   height = 1.0 * 100,
#   units = "mm",
#   dpi = 320
# )
# 
# # Display plot in RStudio
# ric_plot
# 
# 
# # Define base path
# base_dir <- "~/Library/CloudStorage/OneDrive-VrijeUniversiteitAmsterdam/Documents/GitHub/ShrinkageTrees/revision/ric"
# 
# 
# # Save the plot
# ggsave(
#   filename = file.path(base_dir, "revision_ric.png"),
#   plot = ric_plot,
#   width = 1.0*250, height = 0.8*100, units = "mm", dpi = 320
# )
