## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse  = TRUE,
  comment   = "#>",
  message   = FALSE,
  warning   = FALSE,
  fig.width = 6,
  fig.height = 3.5,
  out.width = "100%"
)
library(ShrinkageTrees)
set.seed(42)

## ----tc-data------------------------------------------------------------------
set.seed(50)
n_tc <- 60;  p_tc <- 5
X_tc <- matrix(rnorm(n_tc * p_tc), n_tc, p_tc)
W_tc <- rbinom(n_tc, 1, 0.5)
tau_tc <- 1.5 * (X_tc[, 1] > 0)
y_tc <- X_tc[, 1] + W_tc * tau_tc + rnorm(n_tc, sd = 0.5)

## ----tc-centered--------------------------------------------------------------
# Centered (default)
fit_tc_cen <- CausalHorseForest(
  y = y_tc,
  X_train_control = X_tc, X_train_treat = X_tc,
  treatment_indicator_train = W_tc,
  treatment_coding = "centered",
  number_of_trees = 5, N_post = 50, N_burn = 25,
  store_posterior_sample = TRUE, verbose = FALSE
)
cat("Centered — ATE:",
    round(mean(fit_tc_cen$train_predictions_treat), 3), "\n")

## ----tc-binary----------------------------------------------------------------
# Binary
fit_tc_bin <- CausalHorseForest(
  y = y_tc,
  X_train_control = X_tc, X_train_treat = X_tc,
  treatment_indicator_train = W_tc,
  treatment_coding = "binary",
  number_of_trees = 5, N_post = 50, N_burn = 25,
  store_posterior_sample = TRUE, verbose = FALSE
)
cat("Binary — ATE:",
    round(mean(fit_tc_bin$train_predictions_treat), 3), "\n")

## ----tc-adaptive--------------------------------------------------------------
# Adaptive (requires propensity scores)
ps_tc <- pnorm(0.3 * X_tc[, 1])   # simple propensity model for illustration

fit_tc_ada <- CausalHorseForest(
  y = y_tc,
  X_train_control = X_tc, X_train_treat = X_tc,
  treatment_indicator_train = W_tc,
  treatment_coding = "adaptive",
  propensity = ps_tc,
  number_of_trees = 5, N_post = 50, N_burn = 25,
  store_posterior_sample = TRUE, verbose = FALSE
)
cat("Adaptive — ATE:",
    round(mean(fit_tc_ada$train_predictions_treat), 3), "\n")

## ----tc-invariant-------------------------------------------------------------
# Invariant (parameter-expanded)
fit_tc_inv <- CausalHorseForest(
  y = y_tc,
  X_train_control = X_tc, X_train_treat = X_tc,
  treatment_indicator_train = W_tc,
  treatment_coding = "invariant",
  number_of_trees = 5, N_post = 50, N_burn = 25,
  store_posterior_sample = TRUE, verbose = FALSE
)
cat("Invariant — ATE:",
    round(mean(fit_tc_inv$train_predictions_treat), 3), "\n")

# Posterior draws of b0 and b1 are stored in the fitted object
cat("b0 posterior mean:", round(mean(fit_tc_inv$b0), 3), "\n")
cat("b1 posterior mean:", round(mean(fit_tc_inv$b1), 3), "\n")

## ----load-data----------------------------------------------------------------
library(ShrinkageTrees)
data("pdac")

# Dimensions and column overview
cat("Patients:", nrow(pdac), "\n")
cat("Columns :", ncol(pdac), "\n")
cat("Clinical columns:", paste(names(pdac)[1:13], collapse = ", "), "\n")
cat("Survival: time (months), censoring rate =",
    round(1 - mean(pdac$status), 2), "\n")
cat("Treatment: radiation =", sum(pdac$treatment),
    "/ control =", sum(1 - pdac$treatment), "\n")

## ----prepare-data-------------------------------------------------------------
time      <- pdac$time
status    <- pdac$status
treatment <- pdac$treatment
X         <- as.matrix(pdac[, !(names(pdac) %in% c("time", "status", "treatment"))])

## ----horsetrees-binary, eval=FALSE--------------------------------------------
# ps_fit <- HorseTrees(
#   y            = treatment,
#   X_train      = X,
#   outcome_type = "binary",
#   k            = 0.1,
#   N_post       = 5000,
#   N_burn       = 5000,
#   verbose      = FALSE
# )
# 
# propensity <- pnorm(ps_fit$train_predictions)

## ----horsetrees-binary-small--------------------------------------------------
set.seed(1)
n <- 80;  p <- 10
X_syn  <- matrix(rnorm(n * p), n, p)
W_syn  <- rbinom(n, 1, pnorm(0.8 * X_syn[, 1]))

ps_fit <- HorseTrees(
  y            = W_syn,
  X_train      = X_syn,
  outcome_type = "binary",
  number_of_trees = 5,
  k            = 0.5,
  N_post       = 50,
  N_burn       = 25,
  verbose      = FALSE
)

propensity_syn <- pnorm(ps_fit$train_predictions)
cat("Propensity scores — range: [",
    round(range(propensity_syn), 3), "]\n")

## ----horsetrees-survival------------------------------------------------------
set.seed(2)
log_T <- X_syn[, 1] + rnorm(n)
C     <- rexp(n, 0.5)
y_syn   <- pmin(exp(log_T), C)
d_syn   <- as.integer(exp(log_T) <= C)

ht_surv <- HorseTrees(
  y               = y_syn,
  status          = d_syn,
  X_train         = X_syn,
  outcome_type    = "right-censored",
  timescale       = "time",
  number_of_trees = 5,
  N_post          = 50,
  N_burn          = 25,
  store_posterior_sample = TRUE,
  verbose         = FALSE
)

cat("Posterior mean log-time (first 5 obs):",
    round(ht_surv$train_predictions[1:5], 3), "\n")
cat("Posterior sigma — mean:",
    round(mean(ht_surv$sigma), 3), "\n")

## ----horsetrees-ic------------------------------------------------------------
set.seed(20)

# Generate true event times
true_T <- rexp(n, rate = exp(-0.5 * X_syn[, 1]))

# Create interval-censored observations
left_syn  <- true_T * runif(n, 0.5, 1.0)
right_syn <- true_T * runif(n, 1.0, 1.5)

# Mark some as exact observations and some as right-censored
exact_idx <- sample(n, 25)
left_syn[exact_idx]  <- true_T[exact_idx]
right_syn[exact_idx] <- true_T[exact_idx]

rc_idx <- sample(setdiff(seq_len(n), exact_idx), 15)
right_syn[rc_idx] <- Inf

cat("Exact events:", sum(left_syn == right_syn), "\n")
cat("Interval-censored:", sum(left_syn < right_syn & is.finite(right_syn)), "\n")
cat("Right-censored:", sum(!is.finite(right_syn)), "\n")

ht_ic <- HorseTrees(
  left_time       = left_syn,
  right_time      = right_syn,
  X_train         = X_syn,
  outcome_type    = "interval-censored",
  timescale       = "time",
  number_of_trees = 5,
  N_post          = 50,
  N_burn          = 25,
  store_posterior_sample = TRUE,
  verbose         = FALSE
)

cat("Posterior mean log-time (first 5 obs):",
    round(ht_ic$train_predictions[1:5], 3), "\n")
cat("Posterior sigma — mean:",
    round(mean(ht_ic$sigma), 3), "\n")

## ----sbart-ic-----------------------------------------------------------------
set.seed(21)
fit_sbart_ic <- SurvivalBART(
  left_time       = left_syn,
  right_time      = right_syn,
  X_train         = X_syn,
  number_of_trees = 5,
  N_post          = 50,
  N_burn          = 25,
  verbose         = FALSE
)

cat("SurvivalBART (IC) class:", class(fit_sbart_ic), "\n")

## ----shrinkage-continuous-----------------------------------------------------
set.seed(3)
y_cont <- X_syn[, 1] + 0.5 * X_syn[, 2] + rnorm(n)

# Horseshoe prior (default for HorseTrees)
fit_hs <- ShrinkageTrees(
  y               = y_cont,
  X_train         = X_syn,
  outcome_type    = "continuous",
  prior_type      = "horseshoe",
  local_hp        = 0.1 / sqrt(5),
  global_hp       = 0.1 / sqrt(5),
  number_of_trees = 5,
  N_post          = 50,
  N_burn          = 25,
  verbose         = FALSE
)

# Forest-wide horseshoe (horseshoe_fw)
fit_fw <- ShrinkageTrees(
  y               = y_cont,
  X_train         = X_syn,
  outcome_type    = "continuous",
  prior_type      = "horseshoe_fw",
  local_hp        = 0.1 / sqrt(5),
  global_hp       = 0.1 / sqrt(5),
  number_of_trees = 5,
  N_post          = 50,
  N_burn          = 25,
  verbose         = FALSE
)

cat("Horseshoe   — train RMSE:",
    round(sqrt(mean((fit_hs$train_predictions - y_cont)^2)), 3), "\n")
cat("Horseshoe FW— train RMSE:",
    round(sqrt(mean((fit_fw$train_predictions - y_cont)^2)), 3), "\n")

## ----survival-bart-dart-------------------------------------------------------
set.seed(4)

# SurvivalBART: classical BART prior, AFT likelihood
fit_sbart <- SurvivalBART(
  time            = y_syn,
  status          = d_syn,
  X_train         = X_syn,
  number_of_trees = 5,
  k               = 2.0,
  N_post          = 50,
  N_burn          = 25,
  verbose         = FALSE
)

# SurvivalDART: Dirichlet (DART) splitting prior
fit_sdart <- SurvivalDART(
  time            = y_syn,
  status          = d_syn,
  X_train         = X_syn,
  number_of_trees = 5,
  k               = 2.0,
  N_post          = 50,
  N_burn          = 25,
  verbose         = FALSE
)

cat("SurvivalBART  class:", class(fit_sbart), "\n")
cat("SurvivalDART  class:", class(fit_sdart), "\n")

## ----hd-data------------------------------------------------------------------
set.seed(20)
n_hd <- 60;  p_hd <- 200
X_hd <- matrix(rnorm(n_hd * p_hd), n_hd, p_hd)

# True log-survival depends only on predictors 1, 2, and 3
log_T_hd <- 1.5 * X_hd[, 1] - 1.0 * X_hd[, 2] + 0.5 * X_hd[, 3] + rnorm(n_hd)
C_hd     <- rexp(n_hd, rate = 0.5)
y_hd     <- pmin(exp(log_T_hd), C_hd)
d_hd     <- as.integer(exp(log_T_hd) <= C_hd)

cat("n =", n_hd, "| p =", p_hd,
    "| active predictors = 3",
    "| censoring rate =", round(1 - mean(d_hd), 2), "\n")

## ----hd-horseshoe-------------------------------------------------------------
set.seed(21)
fit_hd_hs <- ShrinkageTrees(
  y               = y_hd,
  status          = d_hd,
  X_train         = X_hd,
  outcome_type    = "right-censored",
  prior_type      = "horseshoe",
  local_hp        = 0.1 / sqrt(10),
  global_hp       = 0.1 / sqrt(10),
  number_of_trees = 10,
  N_post          = 50,
  N_burn          = 25,
  verbose         = FALSE
)

## ----hd-dart------------------------------------------------------------------
set.seed(22)
fit_hd_dart <- SurvivalDART(
  time            = y_hd,
  status          = d_hd,
  X_train         = X_hd,
  number_of_trees = 10,
  rho_dirichlet   = 3,
  N_post          = 50,
  N_burn          = 25,
  verbose         = FALSE
)

## ----hd-compare---------------------------------------------------------------
rmse_hs   <- sqrt(mean((fit_hd_hs$train_predictions  - log_T_hd)^2))
rmse_dart <- sqrt(mean((fit_hd_dart$train_predictions - log_T_hd)^2))

cat(sprintf("%-18s  train RMSE (log-time): %.3f\n", "Horseshoe", rmse_hs))
cat(sprintf("%-18s  train RMSE (log-time): %.3f\n", "DART",      rmse_dart))

## ----hd-vi, eval=requireNamespace("ggplot2", quietly=TRUE)--------------------
plot(fit_hd_dart, type = "vi", n_vi = 10)

## ----causal-data--------------------------------------------------------------
set.seed(5)
tau_true <- 1.5 * (X_syn[, 1] > 0)    # heterogeneous treatment effect
y_causal <- X_syn[, 1] + W_syn * tau_true + rnorm(n, sd = 0.5)

## ----survbcf, eval=FALSE------------------------------------------------------
# # Full analysis (eval=FALSE — use larger MCMC settings in practice)
# fit_sbcf <- SurvivalBCF(
#   time       = time,
#   status     = status,
#   X_train    = X,
#   treatment  = treatment,
#   propensity = propensity,   # from HorseTrees above
#   N_post     = 5000,
#   N_burn     = 5000,
#   verbose    = FALSE
# )

## ----survbcf-small------------------------------------------------------------
set.seed(6)
fit_sbcf <- SurvivalBCF(
  time            = y_syn,
  status          = d_syn,
  X_train         = X_syn,
  treatment       = W_syn,
  number_of_trees_control = 5,
  number_of_trees_treat   = 5,
  N_post          = 50,
  N_burn          = 25,
  verbose         = FALSE
)
cat("SurvivalBCF class:", class(fit_sbcf), "\n")
cat("ATE (posterior mean):",
    round(mean(fit_sbcf$train_predictions_treat), 3), "\n")

## ----survsbcf-small-----------------------------------------------------------
set.seed(7)
fit_ssbcf <- SurvivalShrinkageBCF(
  time            = y_syn,
  status          = d_syn,
  X_train         = X_syn,
  treatment       = W_syn,
  number_of_trees_control = 5,
  number_of_trees_treat   = 5,
  N_post          = 50,
  N_burn          = 25,
  verbose         = FALSE
)
cat("SurvivalShrinkageBCF class:", class(fit_ssbcf), "\n")

## ----causal-horse-------------------------------------------------------------
set.seed(8)
fit_chf <- CausalHorseForest(
  y                         = y_causal,
  X_train_control           = X_syn,
  X_train_treat             = X_syn,
  treatment_indicator_train = W_syn,
  outcome_type              = "continuous",
  number_of_trees           = 5,
  N_post                    = 50,
  N_burn                    = 25,
  store_posterior_sample    = TRUE,
  verbose                   = FALSE
)

cat("CausalHorseForest class:", class(fit_chf), "\n")

# Posterior mean CATE
cate_mean <- fit_chf$train_predictions_treat
cat("CATE — posterior mean (first 5):",
    round(cate_mean[1:5], 3), "\n")

# Posterior ATE
ate_samples <- rowMeans(fit_chf$train_predictions_sample_treat)
cat("ATE posterior mean:",
    round(mean(ate_samples), 3),
    "  95% CI: [",
    round(quantile(ate_samples, 0.025), 3), ",",
    round(quantile(ate_samples, 0.975), 3), "]\n")

## ----causal-horse-test--------------------------------------------------------
set.seed(9)
X_test <- matrix(rnorm(20 * p), 20, p)
W_test <- rbinom(20, 1, 0.5)

fit_chf_test <- CausalHorseForest(
  y                         = y_causal,
  X_train_control           = X_syn,
  X_train_treat             = X_syn,
  treatment_indicator_train = W_syn,
  X_test_control            = X_test,
  X_test_treat              = X_test,
  treatment_indicator_test  = W_test,
  outcome_type              = "continuous",
  number_of_trees           = 5,
  N_post                    = 50,
  N_burn                    = 25,
  store_posterior_sample    = TRUE,
  verbose                   = FALSE
)

cat("Test CATE (first 5):",
    round(fit_chf_test$test_predictions_treat[1:5], 3), "\n")

## ----causal-shrinkage---------------------------------------------------------
set.seed(10)
lh <- 0.1 / sqrt(5)

fit_csf <- CausalShrinkageForest(
  y                         = y_causal,
  X_train_control           = X_syn,
  X_train_treat             = X_syn,
  treatment_indicator_train = W_syn,
  outcome_type              = "continuous",
  prior_type_control        = "horseshoe",
  prior_type_treat          = "horseshoe",
  local_hp_control          = lh,
  global_hp_control         = lh,
  local_hp_treat            = lh,
  global_hp_treat           = lh,
  number_of_trees_control   = 5,
  number_of_trees_treat     = 5,
  N_post                    = 50,
  N_burn                    = 25,
  store_posterior_sample    = TRUE,
  verbose                   = FALSE
)

cat("CausalShrinkageForest class:", class(fit_csf), "\n")
cat("Acceptance ratio (control):",
    round(fit_csf$acceptance_ratio_control, 3), "\n")
cat("Acceptance ratio (treat)  :",
    round(fit_csf$acceptance_ratio_treat, 3), "\n")

## ----causal-shrinkage-fw------------------------------------------------------
set.seed(11)
fit_fw2 <- CausalShrinkageForest(
  y                         = y_causal,
  X_train_control           = X_syn,
  X_train_treat             = X_syn,
  treatment_indicator_train = W_syn,
  outcome_type              = "continuous",
  prior_type_control        = "horseshoe_fw",
  prior_type_treat          = "horseshoe_fw",
  local_hp_control          = lh,
  global_hp_control         = lh,
  local_hp_treat            = lh,
  global_hp_treat           = lh,
  number_of_trees_control   = 5,
  number_of_trees_treat     = 5,
  N_post                    = 50,
  N_burn                    = 25,
  verbose                   = FALSE
)

cat("Forest-wide shrinkage (control, first 5 draws):\n")
print(round(fit_fw2$forestwide_shrinkage_control[1:5], 4))

## ----print--------------------------------------------------------------------
print(fit_chf)

## ----print-csf----------------------------------------------------------------
print(fit_csf)

## ----summary-shrinkage--------------------------------------------------------
smry <- summary(fit_hs)
print(smry)

## ----summary-causal-----------------------------------------------------------
smry_c <- summary(fit_chf)
print(smry_c)

# Access the ATE directly
cat("ATE mean  :", round(smry_c$treatment_effect$ate, 3), "\n")
cat("ATE 95% CI: [",
    round(smry_c$treatment_effect$ate_lower, 3), ",",
    round(smry_c$treatment_effect$ate_upper, 3), "]\n")

## ----summary-pate-mate--------------------------------------------------------
smry_pate <- summary(fit_chf, bayesian_bootstrap = TRUE)   # default
smry_mate <- summary(fit_chf, bayesian_bootstrap = FALSE)

## ----bb-ate-------------------------------------------------------------------
bb <- bayesian_bootstrap_ate(fit_chf)
cat("PATE:", round(bb$pate_mean, 3),
    " 95% CI: [", round(bb$pate_ci$lower, 3), ",",
                  round(bb$pate_ci$upper, 3), "]\n")
cat("MATE:", round(bb$mate_mean, 3),
    " 95% CI: [", round(bb$mate_ci$lower, 3), ",",
                  round(bb$mate_ci$upper, 3), "]\n")

## ----predict------------------------------------------------------------------
X_new  <- matrix(rnorm(10 * p), 10, p)

pred <- predict(fit_hs, newdata = X_new)
print(pred)

## ----predict-ci---------------------------------------------------------------
# Point estimates and 95% credible intervals
head(data.frame(
  mean  = round(pred$mean,  3),
  lower = round(pred$lower, 3),
  upper = round(pred$upper, 3)
))

## ----predict-causal-----------------------------------------------------------
X_new_ctrl  <- matrix(rnorm(10 * p), 10, p)
X_new_treat <- matrix(rnorm(10 * p), 10, p)

pred_c <- predict(fit_chf, newdata_control = X_new_ctrl,
                  newdata_treat = X_new_treat)
print(pred_c)

## ----predict-causal-detail----------------------------------------------------
# Extract individual components
head(data.frame(
  prognostic = round(pred_c$prognostic$mean, 3),
  cate       = round(pred_c$cate$mean, 3),
  total      = round(pred_c$total$mean, 3)
))

## ----plot-trace, eval=requireNamespace("ggplot2", quietly=TRUE)---------------
plot(fit_hs, type = "trace")

## ----plot-ate, eval=requireNamespace("ggplot2", quietly=TRUE)-----------------
plot(fit_chf, type = "ate")                         # PATE (default)
plot(fit_chf, type = "ate", bayesian_bootstrap = FALSE)  # MATE

## ----plot-cate, eval=requireNamespace("ggplot2", quietly=TRUE)----------------
plot(fit_chf, type = "cate")

## ----vi-fit-------------------------------------------------------------------
set.seed(12)
fit_dart <- ShrinkageTrees(
  y               = y_cont,
  X_train         = X_syn,
  outcome_type    = "continuous",
  prior_type      = "dirichlet",
  local_hp        = 0.1 / sqrt(5),
  number_of_trees = 5,
  N_post          = 50,
  N_burn          = 25,
  verbose         = FALSE
)

## ----plot-vi, eval=requireNamespace("ggplot2", quietly=TRUE)------------------
plot(fit_dart, type = "vi", n_vi = 10)

## ----surv-curve-pop, eval=requireNamespace("ggplot2", quietly=TRUE)-----------
# Population-averaged survival curve with 95% credible band
plot(ht_surv, type = "survival")

## ----surv-curve-km, eval=requireNamespace("ggplot2", quietly=TRUE) && requireNamespace("survival", quietly=TRUE)----
# Same curve with the Kaplan-Meier estimate overlaid for comparison
plot(ht_surv, type = "survival", km = TRUE)

## ----surv-curve-ind, eval=requireNamespace("ggplot2", quietly=TRUE)-----------
# Individual survival curves for observations 1, 20, 40, 60, and 80
plot(ht_surv, type = "survival", obs = c(1, 20, 40, 60, 80))

## ----surv-curve-single, eval=requireNamespace("ggplot2", quietly=TRUE)--------
# Single individual with a narrower 90% credible band
plot(ht_surv, type = "survival", obs = 1, level = 0.90)

## ----surv-pred-setup----------------------------------------------------------
# New observations for prediction
set.seed(99)
X_new <- matrix(rnorm(20 * p), ncol = p)
pred_surv <- predict(ht_surv, newdata = X_new)

## ----surv-pred-pop, eval=requireNamespace("ggplot2", quietly=TRUE)------------
# Population-averaged posterior predictive survival curve
plot(pred_surv, type = "survival")

## ----surv-pred-ind, eval=requireNamespace("ggplot2", quietly=TRUE)------------
# Individual posterior predictive curves for selected new observations
plot(pred_surv, type = "survival", obs = c(1, 5, 10))

## ----multi-chain--------------------------------------------------------------
set.seed(13)
fit_2chain <- ShrinkageTrees(
  y               = y_cont,
  X_train         = X_syn,
  outcome_type    = "continuous",
  prior_type      = "horseshoe",
  local_hp        = 0.1 / sqrt(5),
  global_hp       = 0.1 / sqrt(5),
  number_of_trees = 5,
  N_post          = 50,
  N_burn          = 25,
  n_chains        = 2,
  verbose         = FALSE
)

cat("n_chains stored  :", fit_2chain$mcmc$n_chains, "\n")
cat("Total sigma draws:", length(fit_2chain$sigma),
    " (2 chains x 50 draws)\n")
cat("Per-chain acceptance ratios:\n")
print(round(fit_2chain$chains$acceptance_ratios, 3))

## ----multi-chain-causal-------------------------------------------------------
set.seed(14)
fit_causal_2chain <- CausalShrinkageForest(
  y                         = y_causal,
  X_train_control           = X_syn,
  X_train_treat             = X_syn,
  treatment_indicator_train = W_syn,
  outcome_type              = "continuous",
  prior_type_control        = "horseshoe",
  prior_type_treat          = "horseshoe",
  local_hp_control          = lh,
  global_hp_control         = lh,
  local_hp_treat            = lh,
  global_hp_treat           = lh,
  number_of_trees_control   = 5,
  number_of_trees_treat     = 5,
  N_post                    = 50,
  N_burn                    = 25,
  n_chains                  = 2,
  verbose                   = FALSE
)

cat("Pooled sigma draws:", length(fit_causal_2chain$sigma), "\n")
cat("Per-chain acceptance ratios (control):\n")
print(round(fit_causal_2chain$chains$acceptance_ratios_control, 3))

## ----plot-trace-2chain, eval=requireNamespace("ggplot2", quietly=TRUE)--------
plot(fit_2chain, type = "trace")

## ----plot-density-2chain, eval=requireNamespace("ggplot2", quietly=TRUE)------
plot(fit_2chain, type = "density")

## ----coda-summary, eval=requireNamespace("coda", quietly=TRUE)----------------
# summary() includes convergence diagnostics when coda is available
summary(fit_2chain)

## ----coda-diagnostics, eval=requireNamespace("coda", quietly=TRUE)------------
library(coda)
mcmc_obj <- as.mcmc.list(fit_2chain)

# Gelman-Rubin R-hat (values near 1 indicate convergence)
coda::gelman.diag(mcmc_obj)

# Effective sample size
coda::effectiveSize(mcmc_obj)

# Geweke diagnostic (per chain)
coda::geweke.diag(mcmc_obj[[1]])

## ----pdac-ps, eval=FALSE------------------------------------------------------
# data("pdac")
# 
# time      <- pdac$time
# status    <- pdac$status
# treatment <- pdac$treatment
# X         <- as.matrix(pdac[, !(names(pdac) %in% c("time","status","treatment"))])
# 
# set.seed(2025)
# ps_fit <- HorseTrees(
#   y            = treatment,
#   X_train      = X,
#   outcome_type = "binary",
#   k            = 0.1,
#   N_post       = 5000,
#   N_burn       = 5000,
#   verbose      = FALSE
# )
# 
# propensity <- pnorm(ps_fit$train_predictions)

## ----pdac-ps-overlap, eval=FALSE----------------------------------------------
# # Overlap plot
# p0 <- propensity[treatment == 0]
# p1 <- propensity[treatment == 1]
# 
# hist(p0, breaks = 15, col = rgb(1, 0.5, 0, 0.5), xlim = range(propensity),
#      xlab = "Propensity score", main = "Propensity score overlap")
# hist(p1, breaks = 15, col = rgb(0, 0.5, 0, 0.5), add = TRUE)
# legend("topright", legend = c("Control", "Treated"),
#        fill = c(rgb(1,0.5,0,0.5), rgb(0,0.5,0,0.5)))

## ----pdac-causal, eval=FALSE--------------------------------------------------
# # Augment control matrix with propensity scores (BCF-style)
# X_control <- cbind(propensity, X)
# 
# # Log-transform and centre survival times
# log_time <- log(time) - mean(log(time))
# 
# set.seed(2025)
# fit_pdac <- CausalHorseForest(
#   y                         = log_time,
#   status                    = status,
#   X_train_control           = X_control,
#   X_train_treat             = X,
#   treatment_indicator_train = treatment,
#   outcome_type              = "right-censored",
#   timescale                 = "log",
#   number_of_trees           = 200,
#   N_post                    = 5000,
#   N_burn                    = 5000,
#   store_posterior_sample    = TRUE,
#   verbose                   = FALSE
# )

## ----pdac-ate, eval=FALSE-----------------------------------------------------
# # Print model summary
# print(fit_pdac)
# smry_pdac <- summary(fit_pdac)
# print(smry_pdac)
# 
# # ATE
# cat("ATE posterior mean:",
#     round(smry_pdac$treatment_effect$ate, 3), "\n")
# cat("95% CI: [",
#     round(smry_pdac$treatment_effect$ate_lower, 3), ",",
#     round(smry_pdac$treatment_effect$ate_upper, 3), "]\n")

## ----pdac-diag, eval=FALSE----------------------------------------------------
# # Sigma convergence
# plot(fit_pdac, type = "trace")
# 
# # Posterior ATE density
# plot(fit_pdac, type = "ate")
# 
# # CATE caterpillar plot
# plot(fit_pdac, type = "cate")

