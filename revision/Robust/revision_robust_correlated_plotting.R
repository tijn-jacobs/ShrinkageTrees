# Load required libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggh4x)
library(forcats)

# Set base directory to avoid repeated long paths
base_dir <- "~/GitHub/ShrinkageTrees/simulations/Robust sim"

# Load simulation results
combined_results <- readRDS(file.path(base_dir, "revision_robust_correlated_output.rds"))
linear      <- combined_results$results_linear
nonlinear   <- combined_results$results_nonlinear
homo        <- combined_results$results_homo
null        <- combined_results$results_null

# Add SettingType before combining
linear$SettingType    <- "Linear"
nonlinear$SettingType <- "Non-linear"
homo$SettingType      <- "Homogeneous"
null$SettingType      <- "Null"

# Combine all settings
all_results <- bind_rows(linear, nonlinear, homo, null)

# Add squared error for ATE RMSE computation
all_results <- all_results %>%
  mutate(ATE_RMSE = ATE_bias^2)

# List of metrics to summarize
metrics <- c("sigma", "ATE", "ATE_bias", "ATE_RMSE", "ATE_coverage",
             "ATE_CI_Length", "CATE_RPEHE", "CATE_coverage",
             "CATE_CI_Length", "RMSE", "C_Index")

# Summarize metrics by method, dimension, and setting
summary_df <- all_results %>%
  group_by(Method, p_feat, SettingType) %>%
  summarise(across(all_of(metrics), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  mutate(ATE_RMSE = sqrt(ATE_RMSE))

# Clean and pivot for plotting
df_long <- summary_df %>%
  rename(
    ATE_Coverage = ATE_coverage,
    ATE_CILength = ATE_CI_Length,
    CATE_RMSE    = CATE_RPEHE,
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
    SettingType = factor(SettingType, levels = c("Null", "Homogeneous", "Linear", "Non-linear")),
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

# Create the plot
p <- ggplot(df_long, aes(x = p_feat, y = Value, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  geom_hline(
    data = df_long %>% filter(Metric == "Coverage") %>% distinct(Metric),
    aes(yintercept = 0.95),
    linetype = "dashed",
    color = "black"
  ) +
  facet_nested(Metric ~ SettingType + Estimand, scales = "free_y", switch = "y") +
  labs(x = "Number of covariates", y = NULL, fill = "Method", title = NULL) +
  scale_fill_manual(
    values = c("Causal Horseshoe Forest" = "#c7e9c0",
               "AFT-BART-SP" = "#74c476",
               "AFT-BART-NP" = "#238b45")
  ) +
  theme_bw(base_size = 17) +
  theme(
    strip.placement = "outside",
    strip.background = element_rect(fill = "white", color = NA),
    strip.text.x = element_text(face = "bold"),
    text = element_text(family = "Times New Roman"),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Display plot
print(p)

# Save plot
ggsave(
  filename = file.path(base_dir, "robust_correlated_combined_plot.png"),
  plot = p,
  width = 280,
  height = 160,
  units = "mm",
  dpi = 320
)
