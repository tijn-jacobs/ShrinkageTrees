# Load required libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(grid)

# Define base path
base_dir <- "~/Library/CloudStorage/OneDrive-VrijeUniversiteitAmsterdam/Documents/GitHub/ShrinkageTrees/simulations/Deeper Friedman"

# Load simulation result files
low_LD     <- readRDS(file.path(base_dir, "sim_deeper_friedman2_low_LD_output.RDS"))
low_HD     <- readRDS(file.path(base_dir, "sim_deeper_friedman2_low_output.RDS"))
medium_LD  <- readRDS(file.path(base_dir, "sim_deeper_friedman2_medium_LD_output.RDS"))
medium_HD  <- readRDS(file.path(base_dir, "sim_deeper_friedman2_medium_output.RDS"))
high_LD    <- readRDS(file.path(base_dir, "sim_deeper_friedman2_high_LD_output.RDS"))
high_HD    <- readRDS(file.path(base_dir, "sim_deeper_friedman2_high_output.RDS"))

# Combine LD/HD results by signal strength
low    <- rbind(low_LD$results_low, low_HD$results_low)
medium <- rbind(medium_LD$results_medium, medium_HD$results_medium)
high   <- rbind(high_LD$results_high, high_HD$results_high)

# Add ATE_RMSE column (bias^2 -> RMSE)
low    <- low    %>% mutate(ATE_RMSE = ATE_bias^2)
medium <- medium %>% mutate(ATE_RMSE = ATE_bias^2)
high   <- high   %>% mutate(ATE_RMSE = ATE_bias^2)

# Metrics to average
metrics <- c("sigma", "ATE", "ATE_bias", "ATE_RMSE", "ATE_coverage",
             "ATE_CI_Length", "CATE_RPEHE", "CATE_coverage",
             "CATE_CI_Length", "RMSE", "C_Index")

# Compute mean results per method × p
sum_low <- low %>%
  group_by(Method, p_feat) %>%
  summarise(across(all_of(metrics), ~ mean(.x, na.rm = TRUE)),
            .groups = "drop") %>%
  mutate(ATE_RMSE = sqrt(ATE_RMSE))

sum_medium <- medium %>%
  group_by(Method, p_feat) %>%
  summarise(across(all_of(metrics), ~ mean(.x, na.rm = TRUE)),
            .groups = "drop") %>%
  mutate(ATE_RMSE = sqrt(ATE_RMSE))

sum_high <- high %>%
  group_by(Method, p_feat) %>%
  summarise(across(all_of(metrics), ~ mean(.x, na.rm = TRUE)),
            .groups = "drop") %>%
  mutate(ATE_RMSE = sqrt(ATE_RMSE))

# Pivot and annotate each signal strength
long_low <- sum_low %>%
  select(Method, p_feat, ATE_RMSE, CATE_RPEHE) %>%
  pivot_longer(c(ATE_RMSE, CATE_RPEHE), names_to = "Metric",
               values_to = "Value") %>%
  mutate(Metric = recode(Metric, ATE_RMSE = "ATE", CATE_RPEHE = "CATE"),
         Scenario = "Low")

long_medium <- sum_medium %>%
  select(Method, p_feat, ATE_RMSE, CATE_RPEHE) %>%
  pivot_longer(c(ATE_RMSE, CATE_RPEHE), names_to = "Metric",
               values_to = "Value") %>%
  mutate(Metric = recode(Metric, ATE_RMSE = "ATE", CATE_RPEHE = "CATE"),
         Scenario = "Medium")

long_high <- sum_high %>%
  select(Method, p_feat, ATE_RMSE, CATE_RPEHE) %>%
  pivot_longer(c(ATE_RMSE, CATE_RPEHE), names_to = "Metric",
               values_to = "Value") %>%
  mutate(Metric = recode(Metric, ATE_RMSE = "ATE", CATE_RPEHE = "CATE"),
         Scenario = "High")

# Combine all
combined_long <- bind_rows(long_low, long_medium, long_high) %>%
  mutate(
    Scenario = factor(Scenario, levels = c("Low", "Medium", "High")),
    Method = recode(Method,
                    "IndivAFT" = "AFT-BART",
                    "CHF (k=0.1)" = "Causal Horseshoe Forest")
  )

# Plot aesthetics
method_colors <- c("AFT-BART" = "darkorange3",
                   "Causal Horseshoe Forest" = "forestgreen")

linetype_values <- c("ATE" = "dashed", "CATE" = "solid")
shape_values    <- c("ATE" = 1, "CATE" = 0)

# Final plot
final_plot <- ggplot(combined_long,
                     aes(x = p_feat, y = Value,
                         color = Method,
                         linetype = Metric,
                         shape = Metric)) +
  geom_smooth(se = FALSE, linewidth = 0.8) +
  scale_x_log10(breaks = c(10, 100, 1000)) +
  scale_color_manual(values = method_colors) +
  scale_linetype_manual(values = linetype_values) +
  scale_shape_manual(values = shape_values) +
  facet_wrap(~ Scenario, nrow = 1) +
  labs(
    x = "Number of covariates",
    y = "RMSE",
    color = "Method",
    linetype = "Estimand",
    shape = "Estimand"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    text = element_text(family = "Times New Roman"),
    legend.position = "right",
    legend.key.width = unit(1.5, "cm"),
    strip.text.x = element_text(face = "bold"),
    legend.text = element_text(size = 11),
    legend.title = element_text(size = 15),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 15),
    strip.text = element_text(size = 15),
    axis.text.x = element_text(angle = 45, hjust = 1),
  ) +
  guides(
    linetype = guide_legend(override.aes = list(color = "grey40", size = 1.5)),
    shape = guide_legend(override.aes = list(color = "grey40"))
  ) 

print(final_plot)

# Save the plot
ggsave(
  filename = file.path(base_dir, "deeper_friedman.png"),
  plot = final_plot,
  width = 250, height = 100, units = "mm", dpi = 320
)
