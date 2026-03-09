################################################################################
# Generate all figures and outputs for the ShrinkageTrees R Journal paper
#
# Run this script once from the paper/ directory:
#   Rscript generate_figures.R
#
# It saves figures to figures/ and text outputs to outputs/
################################################################################

library(ShrinkageTrees)
library(survival)
library(ggplot2)

dir.create("figures", showWarnings = FALSE)
dir.create("outputs", showWarnings = FALSE)

data("ovarian")

clin      <- ovarian$clinical
time      <- clin$OS_time / 30.44   # convert days to months
status    <- clin$OS_event
treatment <- clin$treatment

X <- cbind(age         = clin$age,
           figo_stage  = clin$figo_stage,
           tumor_grade = clin$tumor_grade,
           ovarian$X)

# --- Data summary ---
sink("outputs/data_summary.txt")
cat("n =", nrow(clin), "| p =", ncol(X),
    "| event rate =", round(mean(status), 2),
    "| carboplatin =", sum(treatment), "/ cisplatin =", sum(1 - treatment), "\n")
sink()

# --- Train / test split ---
set.seed(2025)
train_idx <- sample(seq_len(nrow(clin)), size = floor(0.8 * nrow(clin)))
test_idx  <- setdiff(seq_len(nrow(clin)), train_idx)

X_train <- X[train_idx, ];  X_test <- X[test_idx, ]
time_train <- time[train_idx];  time_test <- time[test_idx]
status_train <- status[train_idx];  status_test <- status[test_idx]

################################################################################
# PART 1: SURVIVAL PREDICTION
################################################################################

cat("=== Fitting HorseTrees ===\n")
fit_horse <- HorseTrees(
  y               = time_train,
  status          = status_train,
  X_train         = X_train,
  X_test          = X_test,
  outcome_type    = "right-censored",
  timescale       = "time",
  number_of_trees = 50,
  k               = 0.1,
  N_post          = 5000,
  N_burn          = 5000,
  n_chains        = 2,
  store_posterior_sample = TRUE,
  verbose         = TRUE
)

# --- Print / summary ---
sink("outputs/print_horse.txt")
print(fit_horse)
sink()

sink("outputs/summary_horse.txt")
summary(fit_horse)
sink()

# --- C-index ---
c_horse <- concordance(Surv(time_test, status_test) ~ fit_horse$test_predictions)

sink("outputs/cindex.txt")
cat("Test C-index — HorseTrees:", round(c_horse$concordance, 3), "\n")
sink()

# --- Convergence diagnostics ---
chain_cols <- c("steelblue", "indianred")

p_trace   <- plot(fit_horse, type = "trace") +
  scale_colour_manual(values = chain_cols)
p_density <- plot(fit_horse, type = "density") +
  scale_colour_manual(values = chain_cols) +
  scale_fill_manual(values = chain_cols)

ggsave("figures/trace_horse.pdf", p_trace, width = 5, height = 3)
ggsave("figures/density_horse.pdf", p_density, width = 5, height = 3)

################################################################################
# PART 2: POSTERIOR SURVIVAL CURVES
################################################################################

pred <- predict(fit_horse, newdata = X_test)

idx <- which.min(abs(pred$mean - median(pred$mean)))

p_individual <- plot(pred, type = "survival", obs = idx)
ggsave("figures/survival_individual.pdf", p_individual, width = 5, height = 3.5)

p_population <- plot(fit_horse, type = "survival", km = TRUE)
ggsave("figures/survival_population.pdf", p_population, width = 5, height = 3.5)

################################################################################
# PART 3: CAUSAL INFERENCE
################################################################################

cat("=== Fitting CausalHorseForest ===\n")

ps_model   <- glm(treatment ~ age + factor(figo_stage) + factor(tumor_grade),
                   family = binomial, data = clin)
propensity <- predict(ps_model, type = "response")
X_control  <- cbind(propensity = propensity, X)

set.seed(2025)
fit_causal <- CausalHorseForest(
  y                         = log(time),
  status                    = status,
  X_train_control           = X_control,
  X_train_treat             = X,
  treatment_indicator_train = treatment,
  outcome_type              = "right-censored",
  timescale                 = "log",
  number_of_trees           = 50,
  k                         = 0.25,
  N_post                    = 5000,
  N_burn                    = 5000,
  n_chains                  = 2,
  store_posterior_sample    = TRUE,
  verbose                   = TRUE
)

sink("outputs/summary_causal.txt")
summary(fit_causal)
sink()

# --- Treatment effect plots ---
p_ate  <- plot(fit_causal, type = "ate")
p_cate <- plot(fit_causal, type = "cate") +
  geom_point(size = 0.8, colour = "steelblue")

ggsave("figures/ate_posterior.pdf", p_ate, width = 4, height = 3)
ggsave("figures/cate_caterpillar.pdf", p_cate, width = 4, height = 3.5)

cat("\n=== All figures and outputs saved ===\n")
cat("Figures: figures/trace_horse.pdf, density_horse.pdf, survival_individual.pdf,\n")
cat("         survival_population.pdf, ate_posterior.pdf, cate_caterpillar.pdf\n")
cat("Outputs: outputs/data_summary.txt, print_horse.txt, summary_horse.txt,\n")
cat("         cindex.txt, patient_profiles.txt, summary_causal.txt\n")




setwd("~/Library/CloudStorage/OneDrive-VrijeUniversiteitAmsterdam/Documents/GitHub/ShrinkageTrees/paper")

rmarkdown::render("motivation-letter/motivation-letter.md")

