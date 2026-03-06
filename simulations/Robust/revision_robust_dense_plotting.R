# Load required libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggh4x)
library(forcats)

# Define base directory
base_dir <- "~/GitHub/ShrinkageTrees/simulations/Robust sim"

# Load simulation results
combined_results <- readRDS(file.path(base_dir, "revision_robust_dense_output.rds"))
homo1 <- combined_results$results_homo1
homo5 <- combined_results$results_homo5
lin1  <- combined_results$results_lin1
lin5  <- combined_results$results_lin5

# Add grouping labels
homo1$SettingType <- "Homogeneous"
homo1$ErrorVar <- "Error variance 1"
homo5$SettingType <- "Homogeneous"
homo5$ErrorVar <- "Error variance 5"
lin1$SettingType  <- "Heterogeneous"
lin1$ErrorVar  <- "Error variance 1"
lin5$SettingType  <- "Heterogeneous"
lin5$ErrorVar  <- "Error variance 5"

# Combine datasets
all_results <- bind_rows(homo1, homo5, lin1, lin5)

# Add squared error for ATE
all_results <- all_results %>%
  mutate(ATE_RMSE = ATE_bias^2)

# Metrics to summarize
metrics <- c("sigma", "ATE", "ATE_bias", "ATE_RMSE", "ATE_coverage",
             "ATE_CI_Length", "CATE_RPEHE", "CATE_coverage",
             "CATE_CI_Length", "RMSE", "C_Index")

# Group and average results
summary_df <- all_results %>%
  group_by(Method, p_feat, SettingType, ErrorVar) %>%
  summarise(across(all_of(metrics), ~ mean(.x, na.rm = TRUE)),
            .groups = "drop") %>%
  mutate(ATE_RMSE = sqrt(ATE_RMSE))

# Prepare long-format data
df_long <- summary_df %>%
  rename(
    ATE_RMSE      = ATE_RMSE,
    ATE_Coverage  = ATE_coverage,
    ATE_CILength  = ATE_CI_Length,
    CATE_RMSE     = CATE_RPEHE,
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
    SettingType = factor(SettingType,
                         levels = c("Homogeneous", "Heterogeneous")),
    ErrorVar = factor(ErrorVar,
                      levels = c("Error variance 1", "Error variance 5")),
    Estimand = factor(Estimand, levels = c("ATE", "CATE")),
    Metric = recode(Metric,
                    "RMSE" = "RMSE",
                    "Coverage" = "Coverage",
                    "CILength" = "Length"),
    Metric = factor(Metric, levels = c("RMSE", "Coverage", "Length")),
    p_feat = as.factor(p_feat),
    Method = factor(Method,
                    levels = c("CHF (k=0.1)", "IndivAFT (SP)", "IndivAFT (NP)")),
    Method = fct_recode(Method,
                        "Causal Horseshoe Forest" = "CHF (k=0.1)",
                        "AFT-BART-SP" = "IndivAFT (SP)",
                        "AFT-BART-NP" = "IndivAFT (NP)")
  )

# Create final plot
p <- ggplot(df_long, aes(x = p_feat, y = Value, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  geom_hline(
    data = df_long %>% filter(Metric == "Coverage") %>%
      distinct(Metric),
    aes(yintercept = 0.95),
    linetype = "dashed", color = "black"
  ) +
  facet_nested(
    Metric ~ SettingType + ErrorVar + Estimand,
    scales = "free_y", switch = "y"
  ) +
  labs(
    x = "Number of covariates", y = NULL, fill = "Method", title = NULL
  ) +
  scale_fill_manual(
    values = c(
      "Causal Horseshoe Forest" = "#c7e9c0",
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

print(p)

# Save plot
ggsave(
  file.path(base_dir, "robust_dense_combined_plot.png"),
  plot = p, width = 280, height = 160, units = "mm", dpi = 320
)
