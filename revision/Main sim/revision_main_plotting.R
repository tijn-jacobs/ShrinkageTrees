# Simulation summary and visualization for linear and non-linear treatment 
# effect scenarios

library(dplyr)
library(ggplot2)
library(ggforce)
library(tidyr)
library(forcats)
library(ggh4x)

# Set base path
base_dir <- "~/Library/CloudStorage/OneDrive-VrijeUniversiteitAmsterdam/Documents/GitHub/ShrinkageTrees/revision/Main sim"

# Load data
linear <- readRDS(file.path(base_dir, "revision_main_linear_low_output.rds"))
nonlinear <- readRDS(file.path(base_dir, "revision_main_nonlinear_low_output.rds"))

results_linear <- linear$results_linear
results_nonlinear <- nonlinear$results_nonlinear

# Add ATE_RMSE as squared bias
add_rmse <- function(df) {
  df %>% mutate(ATE_RMSE = ATE_bias^2)
}

results_linear <- add_rmse(results_linear)
results_nonlinear <- add_rmse(results_nonlinear)

# Define metrics
metrics <- c("sigma", "ATE", "ATE_bias", "ATE_RMSE", "ATE_coverage",
             "ATE_CI_Length", "CATE_RPEHE", "CATE_coverage",
             "CATE_CI_Length", "RMSE", "C_Index")

# Summary function
summarize_df <- function(df) {
  df %>%
    group_by(Method, p_feat) %>%
    summarise(across(all_of(metrics), ~ mean(.x, na.rm = TRUE)),
              .groups = "drop") %>%
    mutate(ATE_RMSE = sqrt(ATE_RMSE))
}

sum_linear <- summarize_df(results_linear)
sum_nonlinear <- summarize_df(results_nonlinear)













library(dplyr)
library(tidyr)
library(ggplot2)
library(ggh4x)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggh4x)


# CV (cross-validated; for non-linear k=0.2?, for linear k=0.5)


# --------------------------------------------------
# Harmonise methods across df_l and df_n
# --------------------------------------------------
df_l <- sum_linear    %>% mutate(Setting = "Linear")
df_n <- sum_nonlinear %>% mutate(Setting = "Non-linear")

df_l <- df_l %>% 
  filter(Method != "CHF_CV")
df_l <- df_l %>%
  mutate(Method = ifelse(Method == "CHF0.5", "CHF_CV", Method))
common_methods <- intersect(unique(df_l$Method), unique(df_n$Method))

df_l <- df_l %>% filter(Method %in% common_methods)
df_n <- df_n %>% filter(Method %in% common_methods)
df   <- bind_rows(df_l, df_n)

#--------------------------------------------------
# 1. Prepare data
#--------------------------------------------------
df_long <- df %>%
  select(Method, p_feat, Setting,
         CATE_RPEHE, CATE_coverage, CATE_CI_Length) %>%
  pivot_longer(
    cols = c(CATE_RPEHE, CATE_coverage, CATE_CI_Length),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Setting = factor(Setting, levels = c("Linear", "Non-linear")),
    Metric = recode(
      Metric,
      CATE_RPEHE     = "RMSE",
      CATE_coverage  = "Coverage",
      CATE_CI_Length = "Length"
    ),
    Metric = factor(Metric, levels = c("RMSE", "Coverage", "Length")),
    p_feat = factor(p_feat, levels = c(100, 1000, 5000)),
    Method = recode(
      Method,
      "CHF_CV"            = "Causal Horseshoe Forest",
      "AFT_BCF"           = "AFT-BCF",
      "AFT_ShrinkageBCF"  = "AFT-Shrinkage BCF",
      "AFT_BART"          = "AFT-BART",
      "AFT_DART"          = "AFT-DART"
    ),
    Method = factor(
      Method,
      levels = c(
        "Causal Horseshoe Forest",
        "AFT-BCF",
        "AFT-Shrinkage BCF",
        "AFT-BART",
        "AFT-DART"
      )
    )
  )

#--------------------------------------------------
# 2. Plot
#--------------------------------------------------
p <- ggplot(df_long, aes(x = p_feat, y = Value, fill = Method)) +
  geom_bar(
    stat = "identity",
    position = position_dodge(width = 0.8),
    color = "black",
    linewidth = 0.6   # <- thicker black outlines
  ) +
  geom_hline(
    data = df_long %>% filter(Metric == "Coverage"),
    aes(yintercept = 0.95),
    linetype = "dashed",
    color = "black",
    linewidth = 0.6
  ) +
  facet_grid2(
    Setting ~ Metric,
    scales = "free_y",
    independent = "y",
    switch = "y"
  ) +
  facetted_pos_scales(
    y = list(
      Metric == "RMSE" ~
        scale_y_continuous(limits = c(0, 1.9),
                           breaks = c(0.5, 1, 1.5),
                           expand = c(0, 0)),
      Metric == "Coverage" ~
        scale_y_continuous(
          limits = c(0, 1.05),
          breaks = c(0.25, 0.5, 0.75, 1),
          expand = c(0, 0)
        ),
      Metric == "Length" ~
        scale_y_continuous(limits = c(0, 4.2),
                           breaks = c(1, 2, 3, 4),
                           expand = c(0, 0)
                           )
    )
  ) +
  labs(
    x = "Number of covariates",
    y = NULL,
    fill = "Method"
  ) +
  scale_fill_manual(values = c(
    "Causal Horseshoe Forest" = "#042E78",  # Deep Blue (Hero)
    "AFT-BCF"                 = "#FFF3E0",  # Very Pale Cream-Orange
    "AFT-Shrinkage BCF"        = "#FFE0B2",  # Pale Peach
    "AFT-BART"                = "#FFCC80",  # Soft Apricot
    "AFT-DART"                = "#FFB74D"   # Mild Orange
  )) + 
  theme_bw(base_size = 21) +
  theme(
    strip.placement   = "outside",
    strip.background  = element_rect(fill = "white", color = NA),
    strip.text.x      = element_text(face = "bold"),
    strip.text.y      = element_text(face = "bold"),
    axis.title.x      = element_text(face = "italic"),
    text              = element_text(family = "Times New Roman"),
    legend.position   = "bottom",
    axis.text.x       = element_text(angle = 45, hjust = 1),
    panel.grid.minor  = element_blank()
  )

print(p)

ggsave(
  filename = file.path(base_dir, "revision_main_sim_plot.png"),
  plot = p,
  width = 1.23 * 250,
  height = 1.23 * 150,
  units = "mm",
  dpi = 320
)
















make_latex_table_cate <- function(df, scenario_label = "linear") {
  
  method_order <- c(
    "CHF_CV",
    "AFT_BCF",
    "AFT_ShrinkageBCF",
    "AFT_BART",
    "AFT_DART"
  )
  
  method_display <- c(
    "CHF_CV"            = "Causal Horseshoe Forest",
    "AFT_BCF"           = "AFT-BCF",
    "AFT_ShrinkageBCF"  = "AFT-ShrinkageBCF",
    "AFT_BART"          = "AFT-BART",
    "AFT_DART"          = "AFT-DART"
  )
  
  p_vals <- sort(unique(df$p_feat))
  
  cat("\\begin{table}[ht]\n")
  cat("\\centering\n")
  cat("\\begin{tabular}{clccc}\n")
  cat("\\toprule\n")
  cat("$p$ & Method & RMSE & Coverage & Length \\\\\n")
  cat("\\midrule\n")
  
  for (p in p_vals) {
    
    sub_df <- df %>%
      filter(p_feat == p) %>%
      filter(Method %in% method_order) %>%
      mutate(Method = factor(Method, levels = method_order)) %>%
      arrange(Method)
    
    for (i in seq_len(nrow(sub_df))) {
      
      row <- sub_df[i, ]
      
      cat(sprintf(
        "%s & %s & %.3f & %.3f & %.3f \\\\\n",
        if (i == 1) p else " ",
        method_display[row$Method],
        row$CATE_RPEHE,
        row$CATE_coverage,
        row$CATE_CI_Length
      ))
    }
    
    cat("\\midrule\n")
  }
  
  cat("\\bottomrule\n")
  cat("\\end{tabular}\n")
  cat(sprintf(
    "\\caption{CATE performance for the %s scenario. Reported are RMSE, coverage, and average 95\\%% credible interval length across methods and dimensions.}\n",
    scenario_label
  ))
  cat(sprintf("\\label{tab:cate_%s}\n", scenario_label))
  cat("\\end{table}\n\n")
}



sum_linear <- sum_linear %>% 
  filter(Method != "CHF_CV")
sum_linear <- sum_linear %>%
  mutate(Method = ifelse(Method == "CHF0.5", "CHF_CV", Method))
make_latex_table_cate(sum_linear,    "linear")
make_latex_table_cate(sum_nonlinear, "nonlinear")
