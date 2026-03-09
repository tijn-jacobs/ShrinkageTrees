################################################################################
# TCGA Ovarian Cancer Analysis with ShrinkageTrees
#
# Data: TCGA-OV (n=357, p=2000 genes + clinical covariates)
#   - OS_time / OS_event: Overall survival (right-censored)
#   - treatment: 1 = carboplatin, 0 = cisplatin
#   - Clinical: age, figo_stage, tumor_grade
#
# Part 1: High-dimensional survival prediction
# Part 2: Causal inference (carboplatin vs cisplatin)
################################################################################

library(ShrinkageTrees)
data(ovarian)

set.seed(2024)

X_gene <- ovarian$X
clin   <- ovarian$clinical
time   <- clin$OS_time
status <- clin$OS_event
n      <- nrow(clin)

cat("=== TCGA Ovarian Cancer Data ===\n")
cat(sprintf("  n = %d, p_gene = %d\n", n, ncol(X_gene)))
cat(sprintf("  Events: %d (%.0f%%), Censored: %d (%.0f%%)\n",
            sum(status), 100 * mean(status),
            sum(1 - status), 100 * mean(1 - status)))
cat(sprintf("  Median OS (uncensored): %.0f days\n", median(time[status == 1])))
cat(sprintf("  Treatment: carboplatin = %d, cisplatin = %d\n",
            sum(clin$treatment == 1), sum(clin$treatment == 0)))

################################################################################
# PART 1: HIGH-DIMENSIONAL SURVIVAL PREDICTION
################################################################################

cat("\n\n")
cat("================================================================\n")
cat("  PART 1: HIGH-DIMENSIONAL SURVIVAL PREDICTION\n")
cat("================================================================\n")

# --- Train / Test split (80/20) ---
train_idx <- sample(seq_len(n), size = floor(0.8 * n))
test_idx  <- setdiff(seq_len(n), train_idx)

X_train <- X_gene[train_idx, ]
X_test  <- X_gene[test_idx, ]

time_train   <- time[train_idx]
status_train <- status[train_idx]
time_test    <- time[test_idx]
status_test  <- status[test_idx]

cat(sprintf("\nTrain: n=%d  Test: n=%d\n", length(train_idx), length(test_idx)))

# --- Evaluation helper ---
compute_cindex <- function(time, status, predicted) {
  if (requireNamespace("survival", quietly = TRUE)) {
    c_obj <- survival::concordance(
      survival::Surv(time, status) ~ predicted
    )
    return(c_obj$concordance)
  }
  return(NA)
}

# MCMC settings
N_post <- 5000
N_burn <- 5000
n_trees <- 50

# --------------------------------------------------------------------------- #
# 1a) SurvivalBART — Classical BART prior (200 trees, standard prior)
# --------------------------------------------------------------------------- #

cat("\n--- 1a) SurvivalBART ---\n")
t0 <- proc.time()

fit_bart <- SurvivalBART(
  time   = time_train,
  status = status_train,
  X_train = X_train,
  X_test  = X_test,
  number_of_trees = n_trees,
  N_post = N_post,
  N_burn = N_burn,
  store_posterior_sample = TRUE,
  verbose = TRUE
)

t_bart <- (proc.time() - t0)[3]

fit_bart
cat("\n")
summary(fit_bart)

cindex_bart_train <- compute_cindex(time_train, status_train, fit_bart$train_predictions)
cindex_bart_test  <- compute_cindex(time_test, status_test, fit_bart$test_predictions)
cat(sprintf("\nSurvivalBART  C-index: train=%.3f  test=%.3f  (%.1fs)\n",
            cindex_bart_train, cindex_bart_test, t_bart))

# --------------------------------------------------------------------------- #
# 1b) SurvivalDART — Dirichlet prior for sparse variable selection
# --------------------------------------------------------------------------- #

cat("\n--- 1b) SurvivalDART ---\n")
t0 <- proc.time()

fit_dart <- SurvivalDART(
  time   = time_train,
  status = status_train,
  X_train = X_train,
  X_test  = X_test,
  number_of_trees = n_trees,
  a_dirichlet = 0.5,
  b_dirichlet = 1.0,
  N_post = N_post,
  N_burn = N_burn,
  store_posterior_sample = TRUE,
  verbose = TRUE
)

t_dart <- (proc.time() - t0)[3]

print(fit_dart)
cat("\n")
print(summary(fit_dart))

cindex_dart_train <- compute_cindex(time_train, status_train, fit_dart$train_predictions)
cindex_dart_test  <- compute_cindex(time_test, status_test, fit_dart$test_predictions)
cat(sprintf("\nSurvivalDART  C-index: train=%.3f  test=%.3f  (%.1fs)\n",
            cindex_dart_train, cindex_dart_test, t_dart))

# --------------------------------------------------------------------------- #
# 1c) HorseTrees — Horseshoe shrinkage prior
# --------------------------------------------------------------------------- #

cat("\n--- 1c) HorseTrees (horseshoe) ---\n")
t0 <- proc.time()

fit_horse <- HorseTrees(
  y      = time_train,
  status = status_train,
  X_train = X_train,
  X_test  = X_test,
  outcome_type = "right-censored",
  timescale = "time",
  number_of_trees = n_trees,
  N_post = N_post,
  N_burn = N_burn,
  store_posterior_sample = TRUE,
  verbose = TRUE
)

t_horse <- (proc.time() - t0)[3]

print(fit_horse)
cat("\n")
print(summary(fit_horse))

cindex_horse_train <- compute_cindex(time_train, status_train, fit_horse$train_predictions)
cindex_horse_test  <- compute_cindex(time_test, status_test, fit_horse$test_predictions)
cat(sprintf("\nHorseTrees    C-index: train=%.3f  test=%.3f  (%.1fs)\n",
            cindex_horse_train, cindex_horse_test, t_horse))

# --------------------------------------------------------------------------- #
# 1d) ShrinkageTrees with Dirichlet + half-Cauchy
# --------------------------------------------------------------------------- #

cat("\n--- 1d) ShrinkageTrees (half-cauchy + dirichlet) ---\n")
t0 <- proc.time()

fit_hc <- ShrinkageTrees(
  y      = time_train,
  status = status_train,
  X_train = X_train,
  X_test  = X_test,
  outcome_type = "right-censored",
  timescale = "time",
  prior_type = "half-cauchy",
  local_hp = 0.1,
  number_of_trees = n_trees,
  N_post = N_post,
  N_burn = N_burn,
  store_posterior_sample = TRUE,
  verbose = TRUE
)

t_hc <- (proc.time() - t0)[3]

print(fit_hc)
cat("\n")
print(summary(fit_hc))

cindex_hc_train <- compute_cindex(time_train, status_train, fit_hc$train_predictions)
cindex_hc_test  <- compute_cindex(time_test, status_test, fit_hc$test_predictions)
cat(sprintf("\nHalf-Cauchy   C-index: train=%.3f  test=%.3f  (%.1fs)\n",
            cindex_hc_train, cindex_hc_test, t_hc))

# --------------------------------------------------------------------------- #
# 1e) Summary comparison table
# --------------------------------------------------------------------------- #

cat("\n\n=== PREDICTION MODEL COMPARISON ===\n")
cat(sprintf("%-20s  %10s  %10s  %8s\n", "Model", "C-train", "C-test", "Time(s)"))
cat(strrep("-", 55), "\n")
cat(sprintf("%-20s  %10.3f  %10.3f  %8.1f\n", "SurvivalBART", cindex_bart_train, cindex_bart_test, t_bart))
cat(sprintf("%-20s  %10.3f  %10.3f  %8.1f\n", "SurvivalDART", cindex_dart_train, cindex_dart_test, t_dart))
cat(sprintf("%-20s  %10.3f  %10.3f  %8.1f\n", "HorseTrees", cindex_horse_train, cindex_horse_test, t_horse))
cat(sprintf("%-20s  %10.3f  %10.3f  %8.1f\n", "Half-Cauchy", cindex_hc_train, cindex_hc_test, t_hc))

# --------------------------------------------------------------------------- #
# 1f) Diagnostic and survival plots (best model)
# --------------------------------------------------------------------------- #

if (requireNamespace("ggplot2", quietly = TRUE)) {
  cat("\n--- Generating plots ---\n")

  # MCMC diagnostics
  p_trace_bart <- plot(fit_bart, type = "trace")
  p_trace_dart <- plot(fit_dart, type = "trace")
  p_trace_horse <- plot(fit_horse, type = "trace")

  # Variable importance (DART only — has Dirichlet prior)
  p_vi_dart <- plot(fit_dart, type = "vi", n_vi = 20)

  # Population-averaged survival curves with Kaplan-Meier overlay
  p_surv_bart  <- plot(fit_bart, type = "survival", km = TRUE)
  p_surv_dart  <- plot(fit_dart, type = "survival", km = TRUE)
  p_surv_horse <- plot(fit_horse, type = "survival", km = TRUE)

  # Individual survival curves for selected test observations
  # Pick 5 observations spanning the predicted risk spectrum
  risk_order <- order(fit_dart$test_predictions)
  obs_selected <- risk_order[round(seq(1, length(risk_order), length.out = 5))]

  pred_dart <- predict(fit_dart, newdata = X_test)
  p_surv_individual <- plot(pred_dart, type = "survival", obs = obs_selected)

  # Save all prediction plots to PDF
  pdf("tests/manual/ovarian_prediction_plots.pdf", width = 10, height = 7)

  print(p_trace_bart + ggplot2::ggtitle("SurvivalBART — Sigma Traceplot"))
  print(p_trace_dart + ggplot2::ggtitle("SurvivalDART — Sigma Traceplot"))
  print(p_trace_horse + ggplot2::ggtitle("HorseTrees — Sigma Traceplot"))

  print(p_vi_dart + ggplot2::ggtitle("SurvivalDART — Top 20 Variable Importance"))

  print(p_surv_bart + ggplot2::ggtitle("SurvivalBART — Population Survival (train)"))
  print(p_surv_dart + ggplot2::ggtitle("SurvivalDART — Population Survival (train)"))
  print(p_surv_horse + ggplot2::ggtitle("HorseTrees — Population Survival (train)"))

  print(p_surv_individual + ggplot2::ggtitle("SurvivalDART — Individual Survival Curves (test)"))

  dev.off()
  cat("Saved: tests/manual/ovarian_prediction_plots.pdf\n")
}


################################################################################
# PART 2: CAUSAL INFERENCE — CARBOPLATIN vs CISPLATIN
################################################################################

cat("\n\n")
cat("================================================================\n")
cat("  PART 2: CAUSAL INFERENCE — CARBOPLATIN vs CISPLATIN\n")
cat("================================================================\n")

treatment <- clin$treatment   # 1 = carboplatin, 0 = cisplatin

cat(sprintf("\n  Carboplatin (W=1): n=%d, events=%d (%.0f%%)\n",
            sum(treatment == 1), sum(status[treatment == 1]),
            100 * mean(status[treatment == 1])))
cat(sprintf("  Cisplatin  (W=0): n=%d, events=%d (%.0f%%)\n",
            sum(treatment == 0), sum(status[treatment == 0]),
            100 * mean(status[treatment == 0])))

# Use all genes + clinical covariates as confounders
X_clinical <- cbind(
  age        = clin$age,
  figo_stage = clin$figo_stage,
  tumor_grade = clin$tumor_grade
)
X_confounders <- cbind(X_clinical, X_gene)

cat(sprintf("  Confounders: %d (3 clinical + %d genes)\n",
            ncol(X_confounders), ncol(X_gene)))

# Propensity scores via logistic regression on clinical covariates
# (use clinical only to avoid overfitting with p >> n)
ps_model <- glm(treatment ~ age + factor(figo_stage) + factor(tumor_grade),
                family = binomial, data = clin)
propensity <- predict(ps_model, type = "response")

cat(sprintf("  Propensity scores: mean=%.3f, range=[%.3f, %.3f]\n",
            mean(propensity), min(propensity), max(propensity)))

# Causal MCMC settings
N_post_causal <- 1000
N_burn_causal <- 500
n_trees_control <- 50
n_trees_treat   <- 25

# --------------------------------------------------------------------------- #
# 2a) SurvivalBCF — Bayesian Causal Forest for survival
# --------------------------------------------------------------------------- #

cat("\n--- 2a) SurvivalBCF (centered coding) ---\n")
t0 <- proc.time()

fit_bcf <- SurvivalBCF(
  time   = time,
  status = status,
  X_train = X_confounders,
  treatment = treatment,
  propensity = propensity,
  treatment_coding = "centered",
  number_of_trees_control = n_trees_control,
  number_of_trees_treat   = n_trees_treat,
  N_post = N_post_causal,
  N_burn = N_burn_causal,
  store_posterior_sample = TRUE,
  verbose = TRUE
)

t_bcf <- (proc.time() - t0)[3]

print(fit_bcf)
cat("\n")
s_bcf <- summary(fit_bcf)
print(s_bcf)

cat(sprintf("\nSurvivalBCF  time: %.1fs\n", t_bcf))
cat(sprintf("  ATE (log-scale): %.4f [%.4f, %.4f]\n",
            s_bcf$treatment_effect$ate,
            s_bcf$treatment_effect$ate_lower,
            s_bcf$treatment_effect$ate_upper))
cat(sprintf("  CATE SD: %.4f\n", s_bcf$treatment_effect$cate_sd))

# --------------------------------------------------------------------------- #
# 2b) SurvivalBCF with adaptive coding (propensity-adjusted)
# --------------------------------------------------------------------------- #

cat("\n--- 2b) SurvivalBCF (adaptive coding) ---\n")
t0 <- proc.time()

fit_bcf_adapt <- SurvivalBCF(
  time   = time,
  status = status,
  X_train = X_confounders,
  treatment = treatment,
  propensity = propensity,
  treatment_coding = "adaptive",
  number_of_trees_control = n_trees_control,
  number_of_trees_treat   = n_trees_treat,
  N_post = N_post_causal,
  N_burn = N_burn_causal,
  store_posterior_sample = TRUE,
  verbose = TRUE
)

t_bcf_adapt <- (proc.time() - t0)[3]

s_bcf_adapt <- summary(fit_bcf_adapt)
cat(sprintf("\nSurvivalBCF (adaptive)  time: %.1fs\n", t_bcf_adapt))
cat(sprintf("  ATE (log-scale): %.4f [%.4f, %.4f]\n",
            s_bcf_adapt$treatment_effect$ate,
            s_bcf_adapt$treatment_effect$ate_lower,
            s_bcf_adapt$treatment_effect$ate_upper))
cat(sprintf("  CATE SD: %.4f\n", s_bcf_adapt$treatment_effect$cate_sd))

# --------------------------------------------------------------------------- #
# 2c) SurvivalShrinkageBCF — Dirichlet sparsity + half-Cauchy shrinkage
# --------------------------------------------------------------------------- #

cat("\n--- 2c) SurvivalShrinkageBCF ---\n")
t0 <- proc.time()

fit_sbcf <- SurvivalShrinkageBCF(
  time   = time,
  status = status,
  X_train = X_confounders,
  treatment = treatment,
  propensity = propensity,
  treatment_coding = "centered",
  a_dir = 0.5,
  b_dir = 1.0,
  number_of_trees_control = n_trees_control,
  number_of_trees_treat   = n_trees_treat,
  N_post = N_post_causal,
  N_burn = N_burn_causal,
  store_posterior_sample = TRUE,
  verbose = TRUE
)

t_sbcf <- (proc.time() - t0)[3]

s_sbcf <- summary(fit_sbcf)
cat(sprintf("\nSurvivalShrinkageBCF  time: %.1fs\n", t_sbcf))
cat(sprintf("  ATE (log-scale): %.4f [%.4f, %.4f]\n",
            s_sbcf$treatment_effect$ate,
            s_sbcf$treatment_effect$ate_lower,
            s_sbcf$treatment_effect$ate_upper))
cat(sprintf("  CATE SD: %.4f\n", s_sbcf$treatment_effect$cate_sd))

# --------------------------------------------------------------------------- #
# 2d) CausalHorseForest — Direct horseshoe causal forest
# --------------------------------------------------------------------------- #

cat("\n--- 2d) CausalHorseForest ---\n")
t0 <- proc.time()

fit_chf <- CausalHorseForest(
  y      = time,
  status = status,
  X_train_control = X_confounders,
  X_train_treat   = X_confounders,
  treatment_indicator_train = treatment,
  outcome_type = "right-censored",
  timescale = "time",
  propensity = propensity,
  treatment_coding = "centered",
  number_of_trees = n_trees_control,
  N_post = N_post_causal,
  N_burn = N_burn_causal,
  store_posterior_sample = TRUE,
  verbose = TRUE
)

t_chf <- (proc.time() - t0)[3]

s_chf <- summary(fit_chf)
cat(sprintf("\nCausalHorseForest  time: %.1fs\n", t_chf))
cat(sprintf("  ATE (log-scale): %.4f [%.4f, %.4f]\n",
            s_chf$treatment_effect$ate,
            s_chf$treatment_effect$ate_lower,
            s_chf$treatment_effect$ate_upper))
cat(sprintf("  CATE SD: %.4f\n", s_chf$treatment_effect$cate_sd))

# --------------------------------------------------------------------------- #
# 2e) Causal model comparison table
# --------------------------------------------------------------------------- #

cat("\n\n=== CAUSAL MODEL COMPARISON (treatment effect on log-survival) ===\n")
cat(sprintf("%-25s  %8s  %8s  %8s  %8s  %8s\n",
            "Model", "ATE", "Lower", "Upper", "CATE_SD", "Time(s)"))
cat(strrep("-", 75), "\n")

models_causal <- list(
  list(name = "SurvivalBCF (centered)",   s = s_bcf,       t = t_bcf),
  list(name = "SurvivalBCF (adaptive)",   s = s_bcf_adapt, t = t_bcf_adapt),
  list(name = "SurvivalShrinkageBCF",     s = s_sbcf,      t = t_sbcf),
  list(name = "CausalHorseForest",        s = s_chf,       t = t_chf)
)

for (m in models_causal) {
  te <- m$s$treatment_effect
  cat(sprintf("%-25s  %8.4f  %8.4f  %8.4f  %8.4f  %8.1f\n",
              m$name, te$ate, te$ate_lower, te$ate_upper, te$cate_sd, m$t))
}

cat("\nInterpretation: ATE > 0 means carboplatin is associated with *longer*\n")
cat("log-survival times (beneficial). The 95% credible interval indicates\n")
cat("statistical significance if it excludes 0.\n")

# --------------------------------------------------------------------------- #
# 2f) Subgroup analysis — CATE by clinical covariates
# --------------------------------------------------------------------------- #

cat("\n\n--- Subgroup Analysis: CATEs by Clinical Characteristics ---\n")

# Use the SurvivalBCF fit (centered) for subgroup analysis
cate <- fit_bcf$train_predictions_treat

# By FIGO stage
cat("\nCATE by FIGO stage (log-scale):\n")
for (stage in sort(unique(clin$figo_stage))) {
  idx <- which(clin$figo_stage == stage)
  cat(sprintf("  Stage %d: mean CATE = %7.4f, sd = %.4f (n=%d)\n",
              stage, mean(cate[idx]), sd(cate[idx]), length(idx)))
}

# By tumor grade
cat("\nCATE by tumor grade (log-scale):\n")
for (grade in sort(unique(clin$tumor_grade))) {
  idx <- which(clin$tumor_grade == grade)
  cat(sprintf("  Grade %d: mean CATE = %7.4f, sd = %.4f (n=%d)\n",
              grade, mean(cate[idx]), sd(cate[idx]), length(idx)))
}

# By age group
age_breaks <- c(0, 55, 65, 75, Inf)
age_labels <- c("<55", "55-64", "65-74", "75+")
age_group  <- cut(clin$age, breaks = age_breaks, labels = age_labels, right = FALSE)

cat("\nCATE by age group (log-scale):\n")
for (lab in age_labels) {
  idx <- which(age_group == lab)
  if (length(idx) > 0) {
    cat(sprintf("  %5s: mean CATE = %7.4f, sd = %.4f (n=%d)\n",
                lab, mean(cate[idx]), sd(cate[idx]), length(idx)))
  }
}

# --------------------------------------------------------------------------- #
# 2g) Causal inference plots
# --------------------------------------------------------------------------- #

if (requireNamespace("ggplot2", quietly = TRUE)) {
  cat("\n--- Generating causal plots ---\n")

  # ATE posterior densities
  p_ate_bcf  <- plot(fit_bcf, type = "ate")
  p_ate_sbcf <- plot(fit_sbcf, type = "ate")
  p_ate_chf  <- plot(fit_chf, type = "ate")

  # CATE caterpillar plots
  p_cate_bcf  <- plot(fit_bcf, type = "cate")
  p_cate_sbcf <- plot(fit_sbcf, type = "cate")

  # Trace plots
  p_trace_bcf <- plot(fit_bcf, type = "trace")

  # Variable importance (ShrinkageBCF has Dirichlet — supports VI)
  p_vi_control <- plot(fit_sbcf, type = "vi", forest = "control", n_vi = 20)
  p_vi_treat   <- plot(fit_sbcf, type = "vi", forest = "treat", n_vi = 20)

  # Save all causal plots to PDF
  pdf("tests/manual/ovarian_causal_plots.pdf", width = 10, height = 7)

  print(p_trace_bcf + ggplot2::ggtitle("SurvivalBCF — Sigma Traceplot"))

  print(p_ate_bcf + ggplot2::ggtitle("SurvivalBCF (centered) — ATE Posterior"))
  print(p_ate_sbcf + ggplot2::ggtitle("SurvivalShrinkageBCF — ATE Posterior"))
  print(p_ate_chf + ggplot2::ggtitle("CausalHorseForest — ATE Posterior"))

  print(p_cate_bcf + ggplot2::ggtitle("SurvivalBCF — CATE Caterpillar Plot"))
  print(p_cate_sbcf + ggplot2::ggtitle("SurvivalShrinkageBCF — CATE Caterpillar Plot"))

  print(p_vi_control + ggplot2::ggtitle("SurvivalShrinkageBCF — Control Forest VI (top 20)"))
  print(p_vi_treat + ggplot2::ggtitle("SurvivalShrinkageBCF — Treatment Forest VI (top 20)"))

  dev.off()
  cat("Saved: tests/manual/ovarian_causal_plots.pdf\n")
}


################################################################################
# PART 3: ADDITIONAL DIAGNOSTICS
################################################################################

cat("\n\n")
cat("================================================================\n")
cat("  PART 3: MCMC DIAGNOSTICS\n")
cat("================================================================\n")

# Multi-chain run for convergence diagnostics (best prediction model)
cat("\n--- Multi-chain SurvivalDART (2 chains) for Gelman-Rubin ---\n")
t0 <- proc.time()

fit_dart_mc <- SurvivalDART(
  time   = time_train,
  status = status_train,
  X_train = X_train,
  X_test  = X_test,
  number_of_trees = n_trees,
  N_post = N_post,
  N_burn = N_burn,
  store_posterior_sample = TRUE,
  n_chains = 2,
  verbose = TRUE
)

t_mc <- (proc.time() - t0)[3]
cat(sprintf("Multi-chain time: %.1fs\n", t_mc))

s_mc <- summary(fit_dart_mc)
print(s_mc)

# coda diagnostics if available
if (requireNamespace("coda", quietly = TRUE)) {
  cat("\n--- coda MCMC Diagnostics ---\n")
  mc_list <- as.mcmc.list(fit_dart_mc)
  cat("Gelman-Rubin diagnostic (sigma):\n")
  print(coda::gelman.diag(mc_list))
  cat("\nEffective sample size (sigma):\n")
  print(coda::effectiveSize(mc_list))
}

# Multi-chain causal for convergence
cat("\n--- Multi-chain SurvivalBCF (2 chains) ---\n")
t0 <- proc.time()

fit_bcf_mc <- SurvivalBCF(
  time   = time,
  status = status,
  X_train = X_confounders,
  treatment = treatment,
  propensity = propensity,
  treatment_coding = "centered",
  number_of_trees_control = n_trees_control,
  number_of_trees_treat   = n_trees_treat,
  N_post = N_post_causal,
  N_burn = N_burn_causal,
  store_posterior_sample = TRUE,
  n_chains = 2,
  verbose = TRUE
)

t_bcf_mc <- (proc.time() - t0)[3]
cat(sprintf("Multi-chain SurvivalBCF time: %.1fs\n", t_bcf_mc))

s_bcf_mc <- summary(fit_bcf_mc)
print(s_bcf_mc)

cat(sprintf("\nATE (multi-chain): %.4f [%.4f, %.4f]\n",
            s_bcf_mc$treatment_effect$ate,
            s_bcf_mc$treatment_effect$ate_lower,
            s_bcf_mc$treatment_effect$ate_upper))

cat("\n\n=== ANALYSIS COMPLETE ===\n")

