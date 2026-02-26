# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)
library(ggh4x)

# Set base directory
base_dir <- "~/GitHub/ShrinkageTrees/simulations/Robust sim"

# Load simulation results
combined_results <- readRDS(file.path(base_dir, "revision_robust_misspecification_output.rds"))
lin_gum  <- combined_results$results_linear_gumbel
non_gum  <- combined_results$results_nonlinear_gumbel
lin_log  <- combined_results$results_linear_logistic
non_log  <- combined_results$results_nonlinear_logistic

# Helper: return cleaned long-format plot data
prepare_plot_data <- function(results_linear, results_nonlinear) {
  metrics <- c("sigma", "ATE", "ATE_bias", "ATE_RMSE", "ATE_coverage",
               "ATE_CI_Length", "CATE_RPEHE", "CATE_coverage",
               "CATE_CI_Length", "RMSE", "C_Index")
  
  # Add squared error
  results_linear <- results_linear %>% mutate(ATE_RMSE = ATE_bias^2)
  results_nonlinear <- results_nonlinear %>% mutate(ATE_RMSE = ATE_bias^2)
  
  # Aggregate
  df_l <- results_linear %>%
    group_by(Method, p_feat) %>%
    summarise(across(all_of(metrics), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
    mutate(ATE_RMSE = sqrt(ATE_RMSE),
           Setting = "Linear")
  
  df_n <- results_nonlinear %>%
    group_by(Method, p_feat) %>%
    summarise(across(all_of(metrics), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
    mutate(ATE_RMSE = sqrt(ATE_RMSE),
           Setting = "Non-linear")
  
  df <- bind_rows(df_l, df_n)
  
  # Clean and pivot
  df_long <- df %>%
    rename(
      ATE_Coverage = ATE_coverage,
      ATE_CILength = ATE_CI_Length,
      CATE_RMSE = CATE_RPEHE,
      CATE_Coverage = CATE_coverage,
      CATE_CILength = CATE_CI_Length
    ) %>%
    pivot_longer(
      cols = c("ATE_RMSE", "ATE_Coverage", "ATE_CILength",
               "CATE_RMSE", "CATE_Coverage", "CATE_CILength"),
      names_to = c("Estimand", "Metric"),
      names_sep = "_",
      values_to = "Value"
    ) %>%
    mutate(
      Setting = factor(Setting, levels = c("Linear", "Non-linear")),
      Estimand = factor(Estimand, levels = c("ATE", "CATE")),
      Metric = recode(Metric,
                      "RMSE" = "RMSE",
                      "Coverage" = "Coverage",
                      "CILength" = "Length"),
      Metric = factor(Metric, levels = c("RMSE", "Coverage", "Length")),
      p_feat = as.factor(p_feat),
      Method = factor(Method, levels = c("CHF (k=0.1)", "IndivAFT (SP)", "IndivAFT (NP)")),
      Method = fct_recode(Method,
                          "Causal Horseshoe Forest" = "CHF (k=0.1)",
                          "AFT-BART-SP" = "IndivAFT (SP)",
                          "AFT-BART-NP" = "IndivAFT (NP)")
    )
  
  return(df_long)
}

# Helper: create plot from long dataframe
make_plot <- function(df_long, title_text) {
  ggplot(df_long, aes(x = p_feat, y = Value, fill = Method)) +
    geom_bar(stat = "identity", position = "dodge", color = "black") +
    geom_hline(
      data = df_long %>% filter(Metric == "Coverage") %>% distinct(Metric),
      aes(yintercept = 0.95),
      linetype = "dashed",
      color = "black"
    ) +
    facet_nested(Metric ~ Setting + Estimand, scales = "free_y", switch = "y") +
    labs(x = "Number of covariates", y = NULL, fill = "Method", title = title_text) +
    scale_fill_manual(
      values = c("Causal Horseshoe Forest" = "#c7e9c0",
                 "AFT-BART-SP" = "#74c476",
                 "AFT-BART-NP" = "#238b45")
    ) +
    theme_bw(base_size = 19) +
    theme(
      strip.placement = "outside",
      strip.background = element_rect(fill = "white", color = NA),
      strip.text.x = element_text(face = "bold"),
      text = element_text(family = "Times New Roman"),
      legend.position = "bottom",
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}

# Plot for logistic 
df_logistic <- prepare_plot_data(lin_log, non_log)
plot_logistic <- make_plot(df_logistic, title_text = NULL)
plot(plot_logistic)
ggsave(
  filename = file.path(base_dir, "robust_logistic_plot.png"),
  plot = plot_logistic,
  width = 1.2*250,
  height = 1.3*150,
  units = "mm",
  dpi = 320
)

# Plot for Gumbel 
df_gumbel <- prepare_plot_data(lin_gum, non_gum)
plot_gumbel <- make_plot(df_gumbel, title_text = NULL)
plot(plot_gumbel)
ggsave(
  filename = file.path(base_dir, "robust_gumbel_plot.png"),
  plot = plot_gumbel,
  width = 1.2*250,
  height = 1.3*150,
  units = "mm",
  dpi = 320
)
