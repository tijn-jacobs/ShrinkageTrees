library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(ggh4x)

# ============================================================
#              METHOD NAMES & COLOR SCHEME
# ============================================================

method_order <- c(
  "CHF", "CHF_CV",
  "AFT_BART", "AFT_BCF", "AFT_DART", "AFT_S_BCF"
)

method_display <- c(
  "CHF"       = "Causal Horseshoe Forest",
  "CHF_CV"    = "Causal Horseshoe Forest (CV)",
  "AFT_BART"  = "AFT-BART",
  "AFT_BCF"   = "AFT-BCF",
  "AFT_DART"  = "AFT-DART",
  "AFT_S_BCF" = "AFT-S-BCF"
)

method_colors <- c(
  "Causal Horseshoe Forest"       = "#6a3d9a",
  "Causal Horseshoe Forest (CV)"  = "#cab2d6",
  "AFT-BART"                      = "#1f78b4",
  "AFT-BCF"                       = "#33a02c",
  "AFT-DART"                      = "#e31a1c",
  "AFT-S-BCF"                     = "#ff7f00"
)

# ============================================================
#              LOAD & SUMMARIZE SIMULATION RESULTS
# ============================================================

base_dir <- "~/Library/CloudStorage/OneDrive-VrijeUniversiteitAmsterdam/Documents/GitHub/ShrinkageTrees/simulations revision"
linear     <- readRDS(file.path(base_dir, "main_linear_output.rds"))
nonlinear  <- readRDS(file.path(base_dir, "main_nonlinear_output.rds"))

results_linear     <- linear$results_linear
results_nonlinear  <- nonlinear$results_nonlinear

# Add RMSE
results_linear     <- results_linear     %>% mutate(ATE_RMSE = ATE_bias^2)
results_nonlinear  <- results_nonlinear  %>% mutate(ATE_RMSE = ATE_bias^2)

metrics <- c("sigma","ATE","ATE_bias","ATE_RMSE","ATE_coverage",
             "ATE_CI_Length","CATE_RPEHE","CATE_coverage",
             "CATE_CI_Length","RMSE","C_Index")

summarize_df <- function(df) {
  df %>%
    group_by(Method, p_feat) %>%
    summarise(across(all_of(metrics), ~ mean(.x, na.rm = TRUE)),
              .groups = "drop") %>%
    mutate(ATE_RMSE = sqrt(ATE_RMSE))
}

sum_linear    <- summarize_df(results_linear)
sum_nonlinear <- summarize_df(results_nonlinear)

df <- bind_rows(
  sum_linear    %>% mutate(Setting = "Linear"),
  sum_nonlinear %>% mutate(Setting = "Non-linear")
)

# ============================================================
#              CLEAN & RESHAPE FOR PLOTTING
# ============================================================

df_long <- df %>%
  mutate(Method = factor(Method, levels = method_order),
         Method = method_display[Method]) %>%
  rename(
    ATE_RMSE      = ATE_RMSE,
    ATE_Coverage  = ATE_coverage,
    ATE_Length    = ATE_CI_Length,
    CATE_RMSE     = CATE_RPEHE,
    CATE_Coverage = CATE_coverage,
    CATE_Length   = CATE_CI_Length
  ) %>%
  pivot_longer(
    cols = c("ATE_RMSE","ATE_Coverage","ATE_Length",
             "CATE_RMSE","CATE_Coverage","CATE_Length"),
    names_to = c("Estimand","Metric"),
    names_sep = "_",
    values_to = "Value"
  ) %>%
  mutate(
    Estimand = factor(Estimand, levels = c("ATE","CATE")),
    Metric   = factor(Metric, levels = c("RMSE","Coverage","Length")),
    Setting  = factor(Setting, levels = c("Linear","Non-linear"))
  )

# ============================================================
#                       FINAL PLOT
# ============================================================

p <- ggplot(df_long, aes(x = factor(p_feat), y = Value, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  geom_hline(
    data = df_long %>% filter(Metric == "Coverage") %>% distinct(Metric),
    aes(yintercept = 0.95),
    linetype = "dashed",
    color = "black"
  ) +
  ggh4x::facet_nested(
    Metric ~ Setting + Estimand,
    scales = "free_y",
    switch = "y"
  ) +
  scale_fill_manual(values = method_colors) +
  labs(
    x = "Number of covariates",
    y = NULL,
    fill = "Method"
  ) +
  theme_bw(base_size = 18) +
  theme(
    strip.placement = "outside",
    strip.background = element_rect(fill = "white", color = NA),
    strip.text.x = element_text(face = "bold"),
    text = element_text(family = "Times New Roman"),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(p)

# Save
ggsave(
  filename = file.path(base_dir, "main_sim_plot_updated.png"),
  plot = p,
  width = 1.4 * 250,
  height = 1.4 * 150,
  units = "mm",
  dpi = 320
)















library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(ggh4x)

# ============================================================
#              UPDATED METHODS (remove AFT-BART, AFT-DART)
# ============================================================

method_order <- c(
  "CHF", "CHF_CV",
  "AFT_BCF", "AFT_S_BCF"
)

method_display <- c(
  "CHF"       = "Causal Horseshoe Forest",
  "CHF_CV"    = "Causal Horseshoe Forest (CV)",
  "AFT_BCF"   = "AFT-BCF",
  "AFT_S_BCF" = "AFT-S-BCF"
)

method_colors <- c(
  "Causal Horseshoe Forest"       = "#6a3d9a",
  "Causal Horseshoe Forest (CV)"  = "#cab2d6",
  "AFT-BCF"                       = "#33a02c",
  "AFT-S-BCF"                     = "#ff7f00"
)

# ============================================================
#              LOAD & SUMMARIZE SIMULATION RESULTS
# ============================================================

base_dir <- "~/Library/CloudStorage/OneDrive-VrijeUniversiteitAmsterdam/Documents/GitHub/ShrinkageTrees/simulations revision"
linear     <- readRDS(file.path(base_dir, "main_linear_output.rds"))
nonlinear  <- readRDS(file.path(base_dir, "main_nonlinear_output.rds"))

results_linear     <- linear$results_linear
results_nonlinear  <- nonlinear$results_nonlinear

# Remove unwanted methods
results_linear     <- results_linear     %>% filter(!(Method %in% c("AFT_BART", "AFT_DART")))
results_nonlinear  <- results_nonlinear  %>% filter(!(Method %in% c("AFT_BART", "AFT_DART")))

# Add RMSE
results_linear     <- results_linear     %>% mutate(ATE_RMSE = ATE_bias^2)
results_nonlinear  <- results_nonlinear  %>% mutate(ATE_RMSE = ATE_bias^2)

metrics <- c("sigma","ATE","ATE_bias","ATE_RMSE","ATE_coverage",
             "ATE_CI_Length","CATE_RPEHE","CATE_coverage",
             "CATE_CI_Length","RMSE","C_Index")

summarize_df <- function(df) {
  df %>%
    group_by(Method, p_feat) %>%
    summarise(across(all_of(metrics), ~ mean(.x, na.rm = TRUE)),
              .groups = "drop") %>%
    mutate(ATE_RMSE = sqrt(ATE_RMSE))
}

sum_linear    <- summarize_df(results_linear)
sum_nonlinear <- summarize_df(results_nonlinear)

df <- bind_rows(
  sum_linear    %>% mutate(Setting = "Linear"),
  sum_nonlinear %>% mutate(Setting = "Non-linear")
)

# ============================================================
#              CLEAN & RESHAPE FOR PLOTTING
# ============================================================

df_long <- df %>%
  mutate(
    Method = factor(Method, levels = method_order),
    Method = method_display[Method]
  ) %>%
  rename(
    ATE_RMSE      = ATE_RMSE,
    ATE_Coverage  = ATE_coverage,
    ATE_Length    = ATE_CI_Length,
    CATE_RMSE     = CATE_RPEHE,
    CATE_Coverage = CATE_coverage,
    CATE_Length   = CATE_CI_Length
  ) %>%
  pivot_longer(
    cols = c("ATE_RMSE","ATE_Coverage","ATE_Length",
             "CATE_RMSE","CATE_Coverage","CATE_Length"),
    names_to = c("Estimand","Metric"),
    names_sep = "_",
    values_to = "Value"
  ) %>%
  mutate(
    Estimand = factor(Estimand, levels = c("ATE","CATE")),
    Metric   = factor(Metric, levels = c("RMSE","Coverage","Length")),
    Setting  = factor(Setting, levels = c("Linear","Non-linear"))
  )

# ============================================================
#                       FINAL PLOT
# ============================================================

p <- ggplot(df_long, aes(x = factor(p_feat), y = Value, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  geom_hline(
    data = df_long %>% filter(Metric == "Coverage") %>% distinct(Metric),
    aes(yintercept = 0.95),
    linetype = "dashed",
    color = "black"
  ) +
  ggh4x::facet_nested(
    Metric ~ Setting + Estimand,
    scales = "free_y",
    switch = "y"
  ) +
  scale_fill_manual(values = method_colors) +
  labs(
    x = "Number of covariates",
    y = NULL,
    fill = "Method"
  ) +
  theme_bw(base_size = 18) +
  theme(
    strip.placement = "outside",
    strip.background = element_rect(fill = "white", color = NA),
    strip.text.x = element_text(face = "bold"),
    text = element_text(family = "Times New Roman"),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(p)

# Save updated figure
ggsave(
  filename = file.path(base_dir, "main_sim_plot_excluding_BART_DART.png"),
  plot = p,
  width = 1.4 * 250,
  height = 1.3 * 150,
  units = "mm",
  dpi = 320
)
