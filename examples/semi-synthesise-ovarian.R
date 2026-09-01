################################################################################
# Semi-synthesise the ovarian dataset
#
# PURPOSE: Replace the real survival outcomes with simulated outcomes from a
#          known DGP, while keeping the real covariates. This gives us:
#          - Known ground truth for the treatment effect (CATE)
#          - Heterogeneous effects across patients
#          NOTE: this script reads the shipped `ovarian` for its real covariates
#          and rebuilds the simulated parts, so it is self-referential. Genes are
#          identified with grep("^ENSG"); an earlier version excluded a list of
#          clinical names instead, which let right_time/left_time enter the gene
#          block and hence the design matrix.
#
#          - Confounding via year_of_diagnosis: cisplatin
#            patients enrolled earlier, carboplatin later — creating a spurious
#            apparent treatment effect that methods must adjust for
#          - Log-normal AFT model calibrated to real TCGA-OV marginal statistics
#          - Right-censored outcome only (interval-censoring added later)
#
# CLINICAL BACKGROUND:
#   Randomized trials (GOG-158, AGO-OVAR-3, Dutch-Danish) showed carbo ≈ cis
#   in efficacy (HR ≈ 1.0). The TCGA-OV observational data shows a spurious
#   time ratio of 1.42 favoring cisplatin, driven by era confounding (cisplatin
#   used pre-2000, carboplatin after) and indication bias. This DGP mimics
#   that structure: a modest true ATE with strong confounding via
#   year_of_diagnosis.
#
# This script constructs the datasets and plots diagnostics only. It fits no
# models; the worked analysis is in
# simulations/r-journal-paper/ovarian_analysis.R.
#
################################################################################

setwd("~/Library/CloudStorage/OneDrive-VrijeUniversiteitAmsterdam/Documents/GitHub/ShrinkageTrees")
devtools::load_all()

library(survival)
library(ggplot2)

set.seed(666)

# =============================================================================
# 1. Load original data and inspect
# =============================================================================

data("ovarian")

# Current format: flat data frame with clinical columns + gene columns
clinical_cols <- c("OS_time", "OS_event", "treatment", "age",
                   "figo_stage", "tumor_grade")

# Remove rows with NA in any clinical covariate (tumor_grade has NAs for GX)
complete <- complete.cases(ovarian[, clinical_cols])
cat(sprintf("Dropping %d rows with missing clinical covariates\n",
            sum(!complete)))
ovarian <- ovarian[complete, ]

clin    <- ovarian[, clinical_cols]

# Identify genes POSITIVELY. Excluding a list of known clinical names is what
# put outcome-derived columns into the design matrix: this script reads the
# shipped `ovarian`, and once right_time/left_time/year_of_diagnosis.1 were in
# it, setdiff() left them in the gene block. Survival times in days have a far
# larger MAD than log2 TPM expression, so they topped the ranking and were
# selected as "genes". grep("^ENSG") cannot pick up a non-gene column.
gene_cols <- grep("^ENSG", names(ovarian), value = TRUE)
X_genes   <- as.matrix(ovarian[, gene_cols])

# Select top 1000 genes by median absolute deviation
gene_mad <- apply(X_genes, 2, mad)
top_genes <- order(gene_mad, decreasing = TRUE)[1:min(1000, ncol(X_genes))]
X_genes <- X_genes[, top_genes]

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

# Median survival (in months)
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
# 2. Generate year_of_diagnosis (confounder)
# =============================================================================
# year_of_diagnosis drives treatment assignment (cisplatin used pre-~2000,
# carboplatin after) and has a direct prognostic effect (time ratio 1.03/year
# from the real Weibull AFT). It is therefore a confounder, not an instrument.
#
# We simulate it to match the TCGA-OV enrollment window (~1992-2013) with
# treatment assignment conditional on year.

n <- nrow(clin)
A_real <- clin$treatment  # 1 = carboplatin, 0 = cisplatin (original assignment)

# Generate year_of_diagnosis: uniform over the TCGA enrollment window,
# with slight correlation to age (older patients slightly later enrollment)
year_center <- 2003  # midpoint of enrollment
age_centered <- clin$age - median(clin$age)
year_of_diagnosis <- round(
  year_center + runif(n, -10, 10) + 0.05 * age_centered + rnorm(n, 0, 1)
)
year_of_diagnosis <- pmax(1992, pmin(2013, year_of_diagnosis))

# Treatment assignment driven by year (the confounding mechanism)
# P(carboplatin = 1) increases with year — logistic model
# Calibrate so ~67% get carbo overall (matching real data: ~258/383)
year_sc <- (year_of_diagnosis - year_center) / 5  # standardise
logit_carbo <- -0.5 + 1.2 * year_sc  # strong year effect
p_carbo <- plogis(logit_carbo)
A <- rbinom(n, 1, p_carbo)

cat(sprintf("\n=== Generated year_of_diagnosis ===\n"))
cat(sprintf("Range: %d – %d, median: %d\n",
            min(year_of_diagnosis), max(year_of_diagnosis),
            median(year_of_diagnosis)))
cat(sprintf("Carboplatin: %d (%.0f%%), Cisplatin: %d (%.0f%%)\n",
            sum(A == 1), 100 * mean(A), sum(A == 0), 100 * (1 - mean(A))))
cat(sprintf("Mean year | carbo: %.1f, cis: %.1f\n",
            mean(year_of_diagnosis[A == 1]), mean(year_of_diagnosis[A == 0])))

# =============================================================================
# 3. Define the data generating process
# =============================================================================
# Log-normal AFT: log(T) = mu(x) + A * tau(x) + sigma * W
# with Gaussian errors, i.e. a log-normal AFT
#
# Calibrated to real TCGA-OV Weibull AFT coefficients:
#   - Intercept ≈ log(45 months) ≈ 3.81 (pooled KM median)
#   - Age: time ratio 0.99/year → coef ≈ -0.008
#   - Year of diagnosis: time ratio 1.03/year → coef ≈ 0.031
#   - FIGO III/IV vs II: modest negative effect
#   - Weibull scale (sigma): 0.647

# Standardise covariates
X_genes_sc <- scale(X_genes)
age_sc     <- (clin$age - median(clin$age))  # centered, in years
figo_late  <- as.integer(clin$figo_stage >= 3)  # III/IV vs II (collapsed)
grade_hi   <- as.integer(clin$tumor_grade >= 3)  # G3/G4 vs G2 (collapsed)
year_sc_prog <- year_of_diagnosis - year_center  # centered at 2003

# --- Prognostic function mu(x) ---
mu_0 <- 3.81  # intercept: log(45 months), pooled median survival

# Clinical effects (calibrated to real AFT)
mu_age   <- -0.008 * age_sc          # time ratio 0.99/year
mu_year  <-  0.031 * year_sc_prog    # time ratio 1.03/year (direct prognostic effect)
mu_figo  <- -0.15  * figo_late       # stage III/IV worse prognosis
mu_grade <- -0.08  * grade_hi        # high grade slightly worse

# Gene effects: 10 active prognostic genes with varying effect sizes
# The last index is the final column of the MAD ranking, whatever its length:
# a literal 1000 goes out of bounds now that the gene block is genuinely genes
# (997) rather than 997 genes plus three contaminants.
active_prog <- c(15, 5, 10, 20, 50, 100, 200, 500, 750, ncol(X_genes_sc))
stopifnot(max(active_prog) <= ncol(X_genes_sc))
beta_prog   <- c(0.12, -0.09, 0.08, -0.07, 0.05,
                 -0.04, 0.03, -0.03, 0.02, -0.02)
mu_genes <- X_genes_sc[, active_prog] %*% beta_prog

# Interaction: gene 1 x figo_stage (biological: gene modifies stage prognosis)
mu_interaction <- 0.06 * X_genes_sc[, 15] * figo_late

mu <- as.numeric(mu_0 + mu_age + mu_year + mu_figo + mu_grade)
                 # + mu_genes + mu_interaction)

cat(sprintf("\nmu(x): mean = %.2f, sd = %.2f, range = [%.2f, %.2f]\n",
            mean(mu), sd(mu), min(mu), max(mu)))

# --- Treatment effect function tau(x) ---
# TRUE CAUSAL EFFECT: mean ≈ 0 (carbo ≈ cis per RCT evidence)
# But heterogeneous: some patients benefit from one agent over the other.
#
# Heterogeneity drivers:
#   - FIGO stage: advanced disease may respond differently
#   - Gene signature: a few genes modify drug response
#   - Age: modest interaction (younger patients tolerate cisplatin better)
#
# The key point: year_of_diagnosis does NOT enter tau. It affects treatment
# assignment and prognosis, but not the treatment effect itself, which is what
# makes it a confounder rather than an effect modifier.

tau_figo <- -0.10 * figo_late          # late stage: slight cis advantage
tau_age  <- -0.003 * age_sc            # older: slight carbo advantage
tau_gene <- -0.15 * X_genes_sc[, 3] +  # gene 3: strong modifier
             0.08 * X_genes_sc[, 7]     # gene 7: moderate modifier

tau <- as.numeric(-0.2 + tau_figo + tau_age + tau_gene)

# Verify: mean should be near 0, sd ≈ 0.15-0.25
cat(sprintf("\ntau(x): mean = %.4f, sd = %.3f, range = [%.3f, %.3f]\n",
            mean(tau), sd(tau), min(tau), max(tau)))
cat(sprintf("Proportion with tau > 0 (benefit from carboplatin): %.2f\n",
            mean(tau > 0)))
hist(tau, breaks = 30, main = "True CATE distribution (check: centered near 0)")
abline(v = 0, col = "red", lwd = 2)

# --- Generate survival times (log-normal AFT) ---
# Log-normal AFT: log(T) = mu + A * tau + sigma * W, with W ~ N(0, 1).
# Deliberately the same error family the package fits: the worked example is
# a well-specified case.
sigma <- 0.65  # Weibull scale from real AFT (was 0.647)
W     <- rnorm(n)  
log_T <- mu + A * tau + sigma * W
T_true <- exp(log_T)

cat(sprintf("\nTrue survival times (months): median = %.1f, mean = %.1f\n",
            median(T_true), mean(T_true)))
cat(sprintf("Median by treatment — carbo: %.1f, cis: %.1f months\n",
            median(T_true[A == 1]), median(T_true[A == 0])))

# =============================================================================
# 4. Apply right-censoring to match target event rate
# =============================================================================

# Target: ~62% events overall, similar rates across treatment groups
# Use an exponential censoring mechanism, tune the rate

generate_censored <- function(T_true, treatment, target_event_rate = 0.32,
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
  list(time = Y, status = delta, censoring_time = C,
       cens_rate_used = cens_rate, achieved_rate = rate)
}

rc <- generate_censored(T_true, A, target_event_rate = 0.62)
time_obs   <- rc$time
status_obs <- rc$status

cat(sprintf("\n=== Right-censored outcome ===\n"))
cat(sprintf("Overall event rate: %.2f (target: 0.62)\n", mean(status_obs)))
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
# 5. Diagnostic plots
# =============================================================================

par(mfrow = c(2, 3))

# 4a. True CATE distribution (should be centered near 0)
hist(tau, breaks = 30, col = "steelblue", border = "white",
     main = "True CATE distribution (log-scale)",
     xlab = expression(tau(x)))
abline(v = mean(tau), col = "red", lwd = 2, lty = 2)
abline(v = 0, col = "black", lwd = 1, lty = 3)
legend("topright", c(paste0("ATE = ", round(mean(tau), 4)), "Zero"),
       col = c("red", "black"), lty = c(2, 3), lwd = c(2, 1))

# 4b. True vs observed survival times
plot(T_true, time_obs, pch = 16, cex = 0.5,
     col = ifelse(status_obs == 1, "black", "red"),
     xlab = "True survival time", ylab = "Observed time",
     main = "True vs observed (right-censored)")
abline(0, 1, col = "grey", lty = 2)
legend("topleft", c("Event", "Censored"), col = c("black", "red"), pch = 16)

# 4c. KM: original vs synthesised
plot(km_orig, col = c("lightblue", "pink"), lwd = 1, lty = 2,
     xlab = "Time (months)", ylab = "Survival",
     main = "Original (dashed) vs synthesised (solid)")
lines(km_synth, col = c("blue", "red"), lwd = 2)
legend("topright",
       c("Orig cisplatin", "Orig carboplatin",
         "Synth cisplatin", "Synth carboplatin"),
       col = c("lightblue", "pink", "blue", "red"),
       lty = c(2, 2, 1, 1), lwd = c(1, 1, 2, 2), cex = 0.7)

# 4d. Confounding: year_of_diagnosis by treatment
boxplot(year_of_diagnosis ~ A, col = c("lightblue", "pink"),
        names = c("Cisplatin (0)", "Carboplatin (1)"),
        main = "Year of diagnosis by treatment (confounding)",
        ylab = "Year of diagnosis")

# 4e. tau vs gene 3 (strongest heterogeneity driver)
plot(X_genes_sc[, 3], tau, pch = 16, cex = 0.5,
     col = adjustcolor("steelblue", 0.6),
     xlab = "Gene 3 expression (standardised)", ylab = expression(tau(x)),
     main = "True CATE vs gene 3 (main modifier)")
abline(h = 0, col = "grey", lty = 2)

# 4f. P(carboplatin) vs year — the confounding relationship
year_grid <- seq(min(year_of_diagnosis), max(year_of_diagnosis), by = 1)
year_grid_sc <- (year_grid - year_center) / 5
p_grid <- plogis(-0.5 + 1.2 * year_grid_sc)
plot(year_grid, p_grid, type = "l", lwd = 2, col = "darkred",
     xlab = "Year of diagnosis", ylab = "P(carboplatin)",
     main = "Treatment propensity vs year (confounding)")
rug(year_of_diagnosis[A == 1], col = "red", side = 3)
rug(year_of_diagnosis[A == 0], col = "blue", side = 1)
legend("topleft", c("Cisplatin", "Carboplatin"),
       col = c("blue", "red"), lty = 1, lwd = 2, cex = 0.7)

par(mfrow = c(1, 1))

# =============================================================================
# 6. Assemble the semi-synthesised dataset
# =============================================================================

# Convert observed time from months to days (to match original OS_time scale)
OS_time_synth  <- round(time_obs * 30.44)
OS_event_synth <- status_obs

# Build new data frame in the NEW flat format (for saving later)
# NOTE: treatment is the *simulated* assignment (driven by year_of_diagnosis),
# not the original real assignment.
ovarian_synth <- data.frame(
  OS_time           = OS_time_synth,
  OS_event          = OS_event_synth,
  treatment         = A,
  age               = clin$age,
  figo_stage        = clin$figo_stage,
  tumor_grade       = clin$tumor_grade,
  year_of_diagnosis = year_of_diagnosis,
  X_genes
)

# Ground truth, for validation. `tau` IS the CATE on the log-time scale, so no
# separate cate column is needed. sigma and the active variable indices are
# attributes rather than columns, being constants.
ovarian_truth_synth <- data.frame(
  mu             = mu,
  tau            = tau,
  f              = mu + A * tau,
  propensity     = p_carbo,
  log_time       = log_T,
  time           = T_true,
  censoring_time = rc$censoring_time
)
attr(ovarian_truth_synth, "sigma")       <- sigma
attr(ovarian_truth_synth, "active_prog") <- active_prog
attr(ovarian_truth_synth, "active_tau")  <- c(3L, 7L)

# The non-gene columns must be exactly the seven we intend to ship. This is the
# check that would have caught the leaked columns.
non_gene <- setdiff(names(ovarian_synth),
                    grep("^ENSG", names(ovarian_synth), value = TRUE))
stopifnot(
  identical(non_gene, c("OS_time", "OS_event", "treatment", "age",
                        "figo_stage", "tumor_grade", "year_of_diagnosis")),
  !any(c("left_time", "right_time", "year_of_diagnosis.1") %in% names(ovarian_synth)),
  nrow(ovarian_truth_synth) == nrow(ovarian_synth)
)

cat(sprintf("\n=== Synthesised dataset ===\n"))
cat(sprintf("n = %d, p = %d  (7 clinical + %d genes)\n",
            nrow(ovarian_synth), ncol(ovarian_synth), ncol(X_genes)))
cat(sprintf("ovarian_truth: %d x %d\n",
            nrow(ovarian_truth_synth), ncol(ovarian_truth_synth)))
cat(sprintf("Event rate: %.2f (target: 0.62)\n", mean(ovarian_synth$OS_event)))
cat(sprintf("Median OS_time (days): %.0f (original: %.0f)\n",
            median(ovarian_synth$OS_time), median(clin$OS_time)))

# =============================================================================
# 7. Save
# =============================================================================
# Guarded, because this overwrites the datasets shipped with the package.
# Run the script normally to inspect the diagnostics, then re-run with --save:
#
#   Rscript examples/semi-synthesise-ovarian.R --save
#
# Or from an interactive session:
#
#   save_data <- TRUE
#   source("examples/semi-synthesise-ovarian.R")
#
# Afterwards: update R/data-documentation.R if the columns changed, then
# devtools::document() and devtools::install(). Every downstream number in the
# manuscript must be regenerated.

if (!exists("save_data"))
  save_data <- "--save" %in% commandArgs(trailingOnly = TRUE)

if (save_data) {
  if (!dir.exists("data"))
    stop("No data/ directory here. Run this from the package root.")

  ovarian       <- ovarian_synth
  ovarian_truth <- ovarian_truth_synth

  save(ovarian,       file = "data/ovarian.rda",       compress = "xz")
  save(ovarian_truth, file = "data/ovarian_truth.rda", compress = "xz")

  cat("\n=== Saved ===\n")
  cat(sprintf("  data/ovarian.rda        %d x %d\n",
              nrow(ovarian), ncol(ovarian)))
  cat(sprintf("  data/ovarian_truth.rda  %d x %d  (sigma = %.3f)\n",
              nrow(ovarian_truth), ncol(ovarian_truth),
              attr(ovarian_truth, "sigma")))
  cat("\nNext: devtools::document(); devtools::install(); restart R.\n")
} else {
  cat("\nNothing saved. Re-run with --save to overwrite data/*.rda\n")
}
# cat("Review the plots. The key checks:\n")
# cat("  1. True ATE should be near 0 (null causal effect per RCT evidence)\n")
# cat("  2. CATE distribution should show spread (heterogeneity exists)\n")
# cat("  3. Naive treatment comparison should show spurious cis advantage\n")
# cat("     (due to year_of_diagnosis confounding)\n")
# cat("  4. Adjusting for year_of_diagnosis should attenuate the effect\n")
# 
