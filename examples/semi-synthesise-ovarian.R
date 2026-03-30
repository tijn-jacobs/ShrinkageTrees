################################################################################
# Semi-synthesise the ovarian dataset
#
# PURPOSE: Replace the real survival outcomes with simulated outcomes from a
#          known DGP, while keeping the real covariates. This gives us:
#          - Known ground truth for the treatment effect (CATE)
#          - More heterogeneity in treatment effects than the real data
#          - Both right-censored and interval-censored versions
#          - Similar marginal statistics to the original data
#
# NOTE: This file is to be deleted once the final dataset is generated.
#
# WORKFLOW:
#   1. Load original data, inspect marginal statistics
#   2. Define DGP: prognostic function mu(x) + heterogeneous tau(x)
#   3. Generate survival times, apply censoring
#   4. Tune censoring to match original event rates (~52% overall)
#   5. Create interval-censored version
#   6. Plot diagnostics (KM curves, CATE distribution, etc.)
#   7. Fit CausalHorseForest, check caterpillar plot and C-index
#   8. If satisfactory, save the dataset
################################################################################
setwd("~/Library/CloudStorage/OneDrive-VrijeUniversiteitAmsterdam/Documents/GitHub/ShrinkageTrees")
devtools::load_all()

library(survival)
library(ggplot2)

set.seed(206)

# =============================================================================
# 1. Load original data and inspect
# =============================================================================

data("ovarian")

# Current format: list with $clinical and $X — extract into working variables
clin    <- ovarian$clinical
X_genes <- ovarian$X

cat("=== Original data ===\n")
cat(sprintf("n = %d, p_genes = %d\n", nrow(clin), ncol(X_genes)))
cat(sprintf("Event rate: %.2f\n", mean(clin$OS_event)))
cat(sprintf("Carboplatin: %d, Cisplatin: %d\n",
            sum(clin$treatment == 1), sum(clin$treatment == 0)))

# Event rate by treatment group
event_carbo <- mean(clin$OS_event[clin$treatment == 1])
event_cis   <- mean(clin$OS_event[clin$treatment == 0])
cat(sprintf("Event rate — carboplatin: %.2f, cisplatin: %.2f\n",
            event_carbo, event_cis))

# Median survival (uncensored, in months)
time_months <- clin$OS_time / 30.44
cat(sprintf("Median OS (events only): %.1f months\n",
            median(time_months[clin$OS_event == 1])))
cat(sprintf("Median OS (all, KM would differ): %.1f months\n",
            median(time_months)))

# KM plot of original data
km_orig <- survfit(Surv(time_months, clin$OS_event) ~ clin$treatment)
plot(km_orig, col = c("blue", "red"), lwd = 2,
     xlab = "Time (months)", ylab = "Survival probability",
     main = "Original data — KM by treatment")
legend("topright", c("Cisplatin (0)", "Carboplatin (1)"),
       col = c("blue", "red"), lty = 1, lwd = 2)

# =============================================================================
# 1b. Fit AFT model to original data (clinical covariates only)
# =============================================================================
# Log-normal AFT to extract baseline parameters for the DGP

aft_fit <- survreg(
  Surv(time_months, clin$OS_event) ~ treatment + age + figo_stage + tumor_grade,
  data = clin,
  dist = "lognormal"
)

cat("\n=== AFT model (log-normal, clinical covariates) ===\n")
print(summary(aft_fit))

# Extract key parameters
cat("\n--- Extracted parameters ---\n")
cat(sprintf("Intercept (baseline log-survival): %.3f\n", coef(aft_fit)["(Intercept)"]))
cat(sprintf("Treatment effect (carboplatin):    %.3f\n", coef(aft_fit)["treatment"]))
cat(sprintf("Age effect:                        %.3f\n", coef(aft_fit)["age"]))
cat(sprintf("FIGO stage effect:                 %.3f\n", coef(aft_fit)["figo_stage"]))
cat(sprintf("Tumor grade effect:                %.3f\n", coef(aft_fit)["tumor_grade"]))
cat(sprintf("Log-scale (sigma):                 %.3f\n", aft_fit$scale))

# =============================================================================
# 2. Define the data generating process
# =============================================================================

# Keep real covariates
n <- nrow(clin)
A <- clin$treatment  # 1 = carboplatin, 0 = cisplatin

# Standardise the gene expression matrix for stable coefficients
X_genes_sc <- scale(X_genes)

# Clinical covariates (standardised)
age_sc   <- scale(clin$age)
figo_sc  <- scale(clin$figo_stage)
grade_sc <- scale(clin$tumor_grade)

# --- Prognostic function mu(x) ---
# Use clinical vars + a handful of genes for the true signal
# Target: median log-survival around log(40 months) ≈ 3.7
mu_0 <- 3.7  # intercept (≈ 40 months median survival)

# Clinical effects
mu_clinical <- -0.15 * age_sc - 0.25 * figo_sc - 0.10 * grade_sc

# Gene effects: pick 10 "active" genes with varying effect sizes
active_prog <- c(1, 5, 10, 20, 50, 100, 200, 500, 1000, 1500)
beta_prog   <- c(0.20, -0.15, 0.12, -0.10, 0.08,
                 -0.07, 0.06, -0.05, 0.04, -0.03)
mu_genes <- X_genes_sc[, active_prog] %*% beta_prog

# Interaction: gene 1 x figo_stage
mu_interaction <- 0.10 * X_genes_sc[, 1] * figo_sc

mu <- as.numeric(mu_0 + mu_clinical + mu_genes + mu_interaction)

cat(sprintf("\nmu(x): mean = %.2f, sd = %.2f, range = [%.2f, %.2f]\n",
            mean(mu), sd(mu), min(mu), max(mu)))

# --- Treatment effect function tau(x) ---
# AFT ATE = -0.22 (carboplatin slightly worse, not significant)
# But with substantial HETEROGENEITY: some patients benefit, others don't
tau_0 <- coef(aft_fit)["treatment"]  # -0.22

# Heterogeneity driven by a few genes and clinical vars
active_tau <- c(3, 15, 75, 300, 800)
beta_tau   <- c(0.12, -0.10, 0.08, -0.06, 0.05)
tau_genes  <- X_genes_sc[, active_tau] %*% beta_tau

# Interaction: treatment effect varies with age and stage
tau_clinical <- 0.15 * age_sc - 0.1 * figo_sc - 0.4 * grade_sc

tau <- as.numeric(tau_0 + tau_genes + tau_clinical)

cat(sprintf("tau(x): mean = %.3f, sd = %.3f, range = [%.3f, %.3f]\n",
            mean(tau), sd(tau), min(tau), max(tau)))
cat(sprintf("Proportion with tau > 0 (benefit from carboplatin): %.2f\n",
            mean(tau > 0)))

# --- Generate log-survival times ---
sigma_eps <- 1.2  # residual SD on log-scale
epsilon   <- rnorm(n, 0, sigma_eps)
log_T     <- mu + A * tau + epsilon
T_true    <- exp(log_T)

cat(sprintf("\nTrue survival times (months): median = %.1f, mean = %.1f\n",
            median(T_true), mean(T_true)))


# =============================================================================
# 3. Apply right-censoring to match original event rate
# =============================================================================

# Target: ~52% events overall, similar rates across treatment groups
# Use an exponential censoring mechanism, tune the rate

generate_censored <- function(T_true, treatment, target_event_rate = 0.52,
                              cens_rate_init = 0.01, tol = 0.01,
                              max_iter = 100) {
  cens_rate <- cens_rate_init
  for (i in seq_len(max_iter)) {
    C <- rexp(length(T_true), rate = cens_rate)
    Y <- pmin(T_true, C)
    delta <- as.integer(T_true <= C)
    rate <- mean(delta)
    if (abs(rate - target_event_rate) < tol) break
    # Adjust: higher cens_rate -> more censoring -> lower event rate
    if (rate > target_event_rate) {
      cens_rate <- cens_rate * 1.1
    } else {
      cens_rate <- cens_rate * 0.9
    }
  }
  list(time = Y, status = delta, cens_rate_used = cens_rate, achieved_rate = rate)
}

rc <- generate_censored(T_true, A, target_event_rate = 0.52)
time_obs   <- rc$time
status_obs <- rc$status

cat(sprintf("\n=== Right-censored outcome ===\n"))
cat(sprintf("Overall event rate: %.2f (target: 0.52)\n", mean(status_obs)))
cat(sprintf("Event rate — carboplatin: %.2f\n",
            mean(status_obs[A == 1])))
cat(sprintf("Event rate — cisplatin:   %.2f\n",
            mean(status_obs[A == 0])))
cat(sprintf("Median observed time: %.1f months\n", median(time_obs)))

# KM plot of synthesised data
km_synth <- survfit(Surv(time_obs, status_obs) ~ A)
plot(km_synth, col = c("blue", "red"), lwd = 2,
     xlab = "Time (months)", ylab = "Survival probability",
     main = "Synthesised data — KM by treatment (right-censored)")
legend("topright", c("Cisplatin (0)", "Carboplatin (1)"),
       col = c("blue", "red"), lty = 1, lwd = 2)


# =============================================================================
# 4. Create interval-censored version
# =============================================================================

# Generate 3 inspection times per patient, uniformly spaced within follow-up
# This creates a mix of: exact, interval-censored, and right-censored

max_followup <- quantile(T_true, 0.95)  # cap at 95th percentile

generate_interval_censored <- function(T_true, max_followup, n_inspections = 3) {
  n <- length(T_true)
  left_time  <- numeric(n)
  right_time <- numeric(n)

  for (i in seq_len(n)) {
    # Random inspection times for this patient
    insp <- sort(runif(n_inspections, min = 0.5, max = max_followup))

    if (T_true[i] <= insp[1]) {
      # Event before first inspection: treat as right-censored at first inspection
      left_time[i]  <- insp[1]
      right_time[i] <- Inf
    } else if (T_true[i] > insp[n_inspections]) {
      # Event after last inspection: right-censored
      left_time[i]  <- insp[n_inspections]
      right_time[i] <- Inf
    } else {
      # Event between two inspections
      idx <- min(which(insp >= T_true[i]))
      left_time[i]  <- insp[idx - 1]
      right_time[i] <- insp[idx]
    }
  }

  list(left_time = left_time, right_time = right_time)
}

ic <- generate_interval_censored(T_true, max_followup)

# Classify observation types
is_exact     <- ic$left_time == ic$right_time  # won't happen with continuous T
is_right_c   <- is.infinite(ic$right_time)
is_interval  <- !is_exact & !is_right_c
is_left_c    <- ic$left_time == 0 & !is_right_c

cat(sprintf("\n=== Interval-censored outcome ===\n"))
cat(sprintf("Left-censored:     %d (%.1f%%)\n",
            sum(is_left_c), 100 * mean(is_left_c)))
cat(sprintf("Interval-censored: %d (%.1f%%)\n",
            sum(is_interval & !is_left_c), 100 * mean(is_interval & !is_left_c)))
cat(sprintf("Right-censored:    %d (%.1f%%)\n",
            sum(is_right_c), 100 * mean(is_right_c)))


# =============================================================================
# 5. Diagnostic plots
# =============================================================================

par(mfrow = c(2, 2))

# 5a. True CATE distribution
hist(tau, breaks = 30, col = "steelblue", border = "white",
     main = "True CATE distribution (log-scale)",
     xlab = expression(tau(x)))
abline(v = mean(tau), col = "red", lwd = 2, lty = 2)
abline(v = 0, col = "black", lwd = 1, lty = 3)
legend("topright", c(paste0("ATE = ", round(mean(tau), 3)), "Zero"),
       col = c("red", "black"), lty = c(2, 3), lwd = c(2, 1))

# 5b. True vs observed survival times
plot(T_true, time_obs, pch = 16, cex = 0.5,
     col = ifelse(status_obs == 1, "black", "red"),
     xlab = "True survival time", ylab = "Observed time",
     main = "True vs observed (right-censored)")
abline(0, 1, col = "grey", lty = 2)
legend("topleft", c("Event", "Censored"), col = c("black", "red"), pch = 16)

# 5c. KM: original vs synthesised
plot(km_orig, col = c("lightblue", "pink"), lwd = 1, lty = 2,
     xlab = "Time (months)", ylab = "Survival",
     main = "Original (dashed) vs synthesised (solid)")
lines(km_synth, col = c("blue", "red"), lwd = 2)
legend("topright",
       c("Orig cisplatin", "Orig carboplatin",
         "Synth cisplatin", "Synth carboplatin"),
       col = c("lightblue", "pink", "blue", "red"),
       lty = c(2, 2, 1, 1), lwd = c(1, 1, 2, 2), cex = 0.7)

# 5d. Interval widths
interval_widths <- ic$right_time - ic$left_time
interval_widths_finite <- interval_widths[is.finite(interval_widths)]
hist(interval_widths_finite, breaks = 30, col = "darkorange", border = "white",
     main = "Interval widths (finite only)",
     xlab = "Width (months)")

par(mfrow = c(1, 1))


# =============================================================================
# 6. Assemble the semi-synthesised dataset
# =============================================================================

# Convert observed time from months to days (to match original OS_time scale)
OS_time_synth  <- round(time_obs * 30.44)
OS_event_synth <- status_obs

# Store interval-censored times in months (the natural scale for the DGP)
# For the dataset we keep days
left_time_days  <- round(ic$left_time * 30.44)
right_time_days <- ifelse(is.infinite(ic$right_time), Inf,
                          round(ic$right_time * 30.44))

# Build new data frame in the NEW flat format (for saving later)
ovarian_synth <- data.frame(
  OS_time     = OS_time_synth,
  OS_event    = OS_event_synth,
  left_time   = left_time_days,
  right_time  = right_time_days,
  treatment   = clin$treatment,
  age         = clin$age,
  figo_stage  = clin$figo_stage,
  tumor_grade = clin$tumor_grade,
  X_genes
)

cat(sprintf("\n=== Synthesised dataset ===\n"))
cat(sprintf("n = %d, p = %d\n", nrow(ovarian_synth), ncol(ovarian_synth)))
cat(sprintf("Event rate: %.2f (original: 0.52)\n", mean(ovarian_synth$OS_event)))
cat(sprintf("Median OS_time (days): %.0f (original: %.0f)\n",
            median(ovarian_synth$OS_time), median(clin$OS_time)))


# =============================================================================
# 7a. Fit single HorseTrees forest + KM comparison
# =============================================================================

cat("\n=== Fitting single HorseTrees forest ===\n")

time_fit   <- ovarian_synth$OS_time / 30.44  # days to months
status_fit <- ovarian_synth$OS_event
treat_fit  <- ovarian_synth$treatment
X_fit      <- as.matrix(ovarian_synth[, -(1:3)])

fit_single <- HorseTrees(
  y               = log(time_fit),
  X_train         = X_fit,
  outcome_type    = "right-censored",
  status          = status_fit,
  timescale       = "log",
  number_of_trees = 200,
  k               = 0.3,
  N_post          = 10,
  N_burn          = 10,
  verbose         = TRUE
)

# C-index for single forest
c_idx_single <- concordance(Surv(time_fit, status_fit) ~ fit_single$train_predictions)
cat(sprintf("C-index (single forest, train): %.3f\n", c_idx_single$concordance))

# KM: observed vs predicted survival
# Predicted median survival time per patient (on original scale)
pred_median <- exp(fit_single$train_predictions)

# Split into risk groups (tertiles of predicted survival)
risk_group <- cut(pred_median,
                  breaks = quantile(pred_median, c(0, 1/3, 2/3, 1)),
                  labels = c("High risk", "Medium risk", "Low risk"),
                  include.lowest = TRUE)

km_risk <- survfit(Surv(time_fit, status_fit) ~ risk_group)
plot(km_risk, col = c("red", "orange", "forestgreen"), lwd = 2,
     xlab = "Time (months)", ylab = "Survival probability",
     main = "KM by predicted risk group (HorseTrees)")
legend("topright", levels(risk_group),
       col = c("red", "orange", "forestgreen"), lty = 1, lwd = 2)

# KM: synthesised data by treatment (compare to original)
km_synth_treat <- survfit(Surv(time_fit, status_fit) ~ treat_fit)
plot(km_orig, col = c("lightblue", "pink"), lwd = 1, lty = 2,
     xlab = "Time (months)", ylab = "Survival probability",
     main = "KM by treatment: original (dashed) vs synthesised (solid)")
lines(km_synth_treat, col = c("blue", "red"), lwd = 2)
legend("topright",
       c("Orig cisplatin", "Orig carboplatin",
         "Synth cisplatin", "Synth carboplatin"),
       col = c("lightblue", "pink", "blue", "red"),
       lty = c(2, 2, 1, 1), lwd = c(1, 1, 2, 2), cex = 0.7)


# =============================================================================
# 7b. Fit CausalHorseForest and evaluate
# =============================================================================

cat("\n=== Fitting CausalHorseForest ===\n")

# Propensity scores (clinical only)
X_ps <- as.matrix(ovarian_synth[, c("age", "figo_stage", "tumor_grade")])
ps_fit <- HorseTrees(
  y               = treat_fit,
  X_train         = X_ps,
  outcome_type    = "binary",
  number_of_trees = 200,
  k               = 0.1,
  N_post          = 2000,
  N_burn          = 2000,
  verbose         = TRUE
)
propensity <- ps_fit$train_predictions

# Fit causal model
X_control <- cbind(propensity = propensity, X_fit)
X_treat   <- X_fit

fit_causal <- CausalHorseForest(
  y                         = log(time_fit),
  status                    = status_fit,
  X_train_control           = X_control,
  X_train_treat             = X_treat,
  treatment_indicator_train = treat_fit,
  outcome_type              = "right-censored",
  timescale                 = "log",
  number_of_trees           = 200,
  k                         = 0.3,
  N_post                    = 10,
  N_burn                    = 10,
  n_chains                  = 1,
  store_posterior_sample     = TRUE,
  verbose                   = TRUE
)

cat("\n")
summary(fit_causal)

# --- C-index ---
cat("\n=== Concordance index ===\n")
# Use the prognostic predictions (mu) for discrimination
c_idx <- concordance(Surv(time_fit, status_fit) ~ fit_causal$train_predictions)
cat(sprintf("C-index (train): %.3f\n", c_idx$concordance))

# --- ATE comparison ---
s <- summary(fit_causal)
cat(sprintf("\nEstimated ATE: %.4f  95%% CI: [%.4f, %.4f]\n",
            s$treatment_effect$ate,
            s$treatment_effect$ate_lower,
            s$treatment_effect$ate_upper))
cat(sprintf("True ATE:      %.4f\n", mean(tau)))

# --- CATE: estimated vs true ---
cate_est <- fit_causal$train_predictions_treat  # posterior mean CATE
cat(sprintf("\nCATE correlation (estimated vs true): %.3f\n", cor(cate_est, tau)))
cat(sprintf("CATE RMSE: %.4f\n", sqrt(mean((cate_est - tau)^2))))

# --- Plots ---
cat("\n=== Generating diagnostic plots ===\n")

# Caterpillar plot
p_cate <- plot(fit_causal, type = "cate")
print(p_cate + ggtitle("Estimated CATEs — semi-synthesised ovarian data"))

# ATE posterior
p_ate <- plot(fit_causal, type = "ate")
print(p_ate + ggtitle("ATE posterior — semi-synthesised ovarian data"))

# Estimated vs true CATE scatter
plot(tau, cate_est, pch = 16, cex = 0.6,
     col = adjustcolor("steelblue", 0.6),
     xlab = "True CATE", ylab = "Estimated CATE",
     main = "Estimated vs true CATE")
abline(0, 1, col = "red", lty = 2, lwd = 2)
abline(h = 0, v = 0, col = "grey", lty = 3)

# Trace plot
p_trace <- plot(fit_causal, type = "trace")
print(p_trace + ggtitle("Sigma traceplot"))


# =============================================================================
# 8. Save if satisfactory
# =============================================================================

# Uncomment the lines below once the dataset looks good.
# This saves in the NEW flat data frame format (not the old list format).
#
ovarian <- ovarian_synth # Make this flat, i.e. just a df not a list of dfs (HERE)
save(ovarian, file = "data/ovarian.rda", compress = "xz")
cat("\nDataset saved to data/ovarian.rda\n")

# Also save the ground truth for validation (not shipped with the package)
ovarian_truth <- data.frame(
  true_log_T = T_true,
  true_mu    = mu,
  true_tau   = tau,
  left_time_days  = left_time_days,
  right_time_days = right_time_days
)
save(ovarian_truth, file = "data/ovarian_truth.rda", compress = "xz")
cat("Ground truth saved to data/ovarian_truth.rda\n")

cat("\n=== Done ===\n")
cat("Review the plots. If the caterpillar plot shows heterogeneity and\n")
cat("the C-index is reasonable, uncomment the save block at the end.\n")

