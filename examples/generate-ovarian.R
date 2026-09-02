###############################################################################
# Regenerate the semi-synthetic `ovarian` and `ovarian_truth` datasets.
#
#   Rscript examples/generate-ovarian.R              # report only
#   Rscript examples/generate-ovarian.R --validate   # + refit the paper's models
#   Rscript examples/generate-ovarian.R --save       # overwrite data/*.rda
#
# Covariates are real: age, FIGO stage, tumour grade and 997 gene expression
# columns come from TCGA-OV. Outcomes, treatment and year of diagnosis are
# simulated from a known process, so the truth is available for validation.
#
# This script is self-referential: it reads the installed `ovarian` to obtain
# the real covariates, then rebuilds the simulated parts. Genes are therefore
# identified positively, with grep("^ENSG"), and never by excluding a list of
# known clinical names. The previous version excluded a six-name list, so when
# `right_time`, `left_time` and `year_of_diagnosis.1` were added to the data
# they were treated as genes; survival times in days have a far larger MAD than
# log2 TPM expression, so they topped the ranking and were selected into the
# gene block. That is how outcome-derived columns entered the design matrix.
#
# Requires ShrinkageTrees (for the input data only).
###############################################################################

library(ShrinkageTrees)
library(survival)

set.seed(666)

cli       <- commandArgs(trailingOnly = TRUE)
save_data <- "--save" %in% cli
validate  <- "--validate" %in% cli

# ---------------------------------------------------------------------------
# 1. Real covariates
# ---------------------------------------------------------------------------

data("ovarian")

clinical_cols <- c("age", "figo_stage", "tumor_grade")
outcome_cols  <- c("OS_time", "OS_event", "treatment")

# Positive identification: anything that is not an ENSG column is not a gene.
gene_cols <- grep("^ENSG", names(ovarian), value = TRUE)

complete <- complete.cases(ovarian[, c(clinical_cols, outcome_cols)])
cat("Dropping ", sum(!complete), " rows with missing clinical covariates\n",
    sep = "")
src <- ovarian[complete, ]

clin    <- src[, clinical_cols]
X_genes <- as.matrix(src[, gene_cols])

# Rank by median absolute deviation, as before. Only genes are eligible.
gene_mad  <- apply(X_genes, 2, mad)
top_genes <- order(gene_mad, decreasing = TRUE)[seq_len(min(1000, ncol(X_genes)))]
X_genes   <- X_genes[, top_genes]

n <- nrow(clin)
cat("Real covariates: n = ", n, ", genes = ", ncol(X_genes), "\n", sep = "")

# ---------------------------------------------------------------------------
# 2. year_of_diagnosis, the confounder
# ---------------------------------------------------------------------------
# It drives treatment assignment (cisplatin pre-~2000, carboplatin after) and
# has a direct prognostic effect, which is what makes it a confounder rather
# than an instrument.

year_center  <- 2003
age_centered <- clin$age - median(clin$age)
year_of_diagnosis <- round(
  year_center + runif(n, -10, 10) + 0.05 * age_centered + rnorm(n, 0, 1)
)
year_of_diagnosis <- pmax(1992, pmin(2013, year_of_diagnosis))

year_sc     <- (year_of_diagnosis - year_center) / 5
logit_carbo <- -0.5 + 1.2 * year_sc
p_carbo     <- plogis(logit_carbo)
A           <- rbinom(n, 1, p_carbo)

cat("Carboplatin: ", sum(A == 1), " (", round(100 * mean(A)), "%)\n", sep = "")

# ---------------------------------------------------------------------------
# 3. Data-generating process
# ---------------------------------------------------------------------------
# Weibull AFT: log(T) = mu(x) + A * tau(x) + sigma * W
# Calibrated to the real TCGA-OV Weibull AFT.

X_genes_sc   <- scale(X_genes)
age_sc       <- clin$age - median(clin$age)
figo_late    <- as.integer(clin$figo_stage >= 3)
grade_hi     <- as.integer(clin$tumor_grade >= 3)
year_sc_prog <- year_of_diagnosis - year_center

# Prognostic function mu(x)
mu_0     <- 3.81                     # log(45 months), pooled median survival
mu_age   <- -0.008 * age_sc          # time ratio 0.99/year
mu_year  <-  0.031 * year_sc_prog    # time ratio 1.03/year
mu_figo  <- -0.15  * figo_late
mu_grade <- -0.08  * grade_hi

# Gene effects are computed but deliberately NOT added to mu, so the genes are
# noise in the outcome model. Enabling them is a decision for a later release;
# it changes every downstream result.
active_prog <- c(15, 5, 10, 20, 50, 100, 200, 500, 750, 1000)
beta_prog   <- c(0.12, -0.09, 0.08, -0.07, 0.05,
                 -0.04, 0.03, -0.03, 0.02, -0.02)
mu_genes       <- X_genes_sc[, active_prog] %*% beta_prog
mu_interaction <- 0.06 * X_genes_sc[, 15] * figo_late

mu <- as.numeric(mu_0 + mu_age + mu_year + mu_figo + mu_grade)
                 # + mu_genes + mu_interaction)

# Treatment effect tau(x): heterogeneous, mean near zero. year_of_diagnosis
# does not enter, so it confounds without modifying the effect.
tau_figo <- -0.10  * figo_late
tau_age  <- -0.003 * age_sc
tau_gene <- -0.15 * X_genes_sc[, 3] + 0.08 * X_genes_sc[, 7]
tau      <- as.numeric(-0.2 + tau_figo + tau_age + tau_gene)

sigma  <- 0.65
W      <- rnorm(n)
log_T  <- mu + A * tau + sigma * W
T_true <- exp(log_T)

cat("mu:  mean ", round(mean(mu), 2), ", sd ", round(sd(mu), 2), "\n", sep = "")
cat("tau: mean ", round(mean(tau), 4), ", sd ", round(sd(tau), 3), "\n", sep = "")

# ---------------------------------------------------------------------------
# 4. Right-censoring
# ---------------------------------------------------------------------------

generate_censored <- function(T_true, target_event_rate = 0.62,
                              cens_rate_init = 0.01, tol = 0.01,
                              max_iter = 100) {
  cens_rate <- cens_rate_init
  for (i in seq_len(max_iter)) {
    C     <- rexp(length(T_true), rate = cens_rate)
    Y     <- pmin(T_true, C)
    delta <- as.integer(T_true <= C)
    rate  <- mean(delta)
    if (abs(rate - target_event_rate) < tol) break
    cens_rate <- if (rate > target_event_rate) cens_rate * 1.1 else cens_rate * 0.9
  }
  list(time = Y, status = delta, censoring_time = C, achieved_rate = rate)
}

rc <- generate_censored(T_true, target_event_rate = 0.62)
cat("Event rate: ", round(rc$achieved_rate, 3), " (target 0.62)\n", sep = "")

# ---------------------------------------------------------------------------
# 5. Assemble
# ---------------------------------------------------------------------------
# Baseline covariates only. left_time and right_time are deliberately absent:
# they are deterministic in (OS_time, OS_event), so they carry no information
# and can only leak the outcome. Rebuild them where interval censoring is
# demonstrated:
#   left  <- ovarian$OS_time
#   right <- ifelse(ovarian$OS_event == 1, ovarian$OS_time, Inf)

ovarian_new <- data.frame(
  OS_time           = rc$time,
  OS_event          = rc$status,
  treatment         = A,
  age               = clin$age,
  figo_stage        = clin$figo_stage,
  tumor_grade       = clin$tumor_grade,
  year_of_diagnosis = year_of_diagnosis,
  X_genes
)

ovarian_truth_new <- data.frame(
  mu             = mu,
  tau            = tau,
  propensity     = p_carbo,
  log_time       = log_T,
  time           = T_true,
  censoring_time = rc$censoring_time
)
attr(ovarian_truth_new, "sigma")       <- sigma
attr(ovarian_truth_new, "active_prog") <- active_prog
attr(ovarian_truth_new, "active_tau")  <- c(3L, 7L)

# ---------------------------------------------------------------------------
# 6. Checks
# ---------------------------------------------------------------------------

n_gene <- length(grep("^ENSG", names(ovarian_new)))
stopifnot(
  nrow(ovarian_new) == n,
  nrow(ovarian_truth_new) == n,
  ncol(ovarian_new) == 7L + n_gene,
  !any(c("left_time", "right_time", "year_of_diagnosis.1") %in% names(ovarian_new)),
  # every non-gene column is one of the seven we intend to ship
  setdiff(names(ovarian_new), grep("^ENSG", names(ovarian_new), value = TRUE)) ==
    c("OS_time", "OS_event", "treatment", "age", "figo_stage",
      "tumor_grade", "year_of_diagnosis")
)

cat("\novarian:       ", nrow(ovarian_new), " x ", ncol(ovarian_new),
    "  (7 clinical + ", n_gene, " genes)\n", sep = "")
cat("ovarian_truth: ", nrow(ovarian_truth_new), " x ",
    ncol(ovarian_truth_new), "\n", sep = "")
cat("design matrix would be: ", n, " x ", 4L + n_gene, "\n", sep = "")
cat("\nTrue ATE (log scale): ", round(mean(tau), 4), "\n", sep = "")
cat("Median true survival: ", round(median(T_true), 1), " months\n", sep = "")

# ---------------------------------------------------------------------------
# 7. Validation: the paper's analysis, run on the data just generated
# ---------------------------------------------------------------------------
# Same models and same design matrix as Sections 4.1 and 4.2 of the
# manuscript, at N_post = N_burn = 1000 so this stays a check rather than a
# rerun of the paper. Because the truth is known here, the estimates can be
# compared against it directly.

if (validate) {

  time_m <- ovarian_new$OS_time
  status <- ovarian_new$OS_event
  trt    <- ovarian_new$treatment

  drop <- c("OS_time", "OS_event", "treatment")
  X    <- as.matrix(ovarian_new[, setdiff(names(ovarian_new), drop)])
  cat("\nDesign matrix: ", nrow(X), " x ", ncol(X), "\n", sep = "")

  n_post <- 1000L
  n_burn <- 1000L

  # -- Section 4.1: standard BART baseline ---------------------------------
  set.seed(42)
  fit_bart <- SurvivalBART(
    time            = time_m,
    status          = status,
    X_train         = X,
    timescale       = "time",
    number_of_trees = 200,
    N_post          = n_post,
    N_burn          = n_burn,
    n_chains        = 4,
    verbose         = FALSE
  )
  c_bart <- concordance(Surv(time_m, status) ~ fit_bart$train_predictions)

  # -- Section 4.1: Horseshoe Forest ---------------------------------------
  set.seed(42)
  fit_horse <- HorseTrees(
    y               = time_m,
    status          = status,
    X_train         = X,
    outcome_type    = "right-censored",
    timescale       = "time",
    number_of_trees = 200,
    k               = 1,
    N_post          = n_post,
    N_burn          = n_burn,
    n_chains        = 4,
    verbose         = FALSE
  )
  c_horse <- concordance(Surv(time_m, status) ~ fit_horse$train_predictions)

  # The DGP caps concordance: with signal sd s and noise sigma, the best any
  # model can reach is 1/2 + asin(s / sqrt(s^2 + sigma^2)) / pi.
  eta   <- mu + A * tau
  rho   <- sd(eta) / sqrt(var(eta) + sigma^2)
  c_max <- 0.5 + asin(rho) / pi

  cat("\n--- C-index ---\n")
  cat("  SurvivalBART   : ", round(c_bart$concordance, 3), "\n", sep = "")
  cat("  HorseTrees     : ", round(c_horse$concordance, 3), "\n", sep = "")
  cat("  DGP ceiling    : ", round(c_max, 3),
      "  (no model can beat this)\n", sep = "")
  if (max(c_bart$concordance, c_horse$concordance) > c_max + 0.05)
    cat("  WARNING: above the ceiling, which indicates leakage in X.\n")

  # -- Section 4.2: propensity and causal forest ---------------------------
  set.seed(42)
  ps_fit <- HorseTrees(
    y               = trt,
    X_train         = X,
    outcome_type    = "binary",
    number_of_trees = 200,
    k               = 1,
    N_post          = n_post,
    N_burn          = n_burn,
    verbose         = FALSE
  )
  propensity <- ps_fit$train_predictions
  X_control  <- cbind(propensity = propensity, X)

  set.seed(42)
  fit_causal <- CausalHorseForest(
    y                         = log(time_m),
    status                    = status,
    X_train_control           = X_control,
    X_train_treat             = X,
    treatment_indicator_train = trt,
    outcome_type              = "right-censored",
    timescale                 = "log",
    number_of_trees           = 200,
    k                         = 1.5,
    N_post                    = n_post,
    N_burn                    = n_burn,
    store_posterior_sample    = TRUE,
    n_chains                  = 4,
    verbose                   = FALSE
  )
  cs  <- summary(fit_causal)
  te  <- cs$treatment_effect
  cate_hat <- fit_causal$train_predictions_treat

  cat("\n--- Causal estimates against the truth ---\n")
  cat("  true ATE       : ", round(mean(tau), 4), "\n", sep = "")
  cat("  estimated ATE  : ", round(te$ate, 4),
      "  95% CI [", round(te$ate_lower, 3), ", ",
      round(te$ate_upper, 3), "]\n", sep = "")
  cat("  bias           : ", round(te$ate - mean(tau), 4), "\n", sep = "")
  cat("  ATE covered    : ",
      mean(tau) >= te$ate_lower && mean(tau) <= te$ate_upper, "\n", sep = "")
  cat("  CATE RMSE      : ", round(sqrt(mean((cate_hat - tau)^2)), 4),
      "\n", sep = "")
  cat("  CATE corr      : ", round(cor(cate_hat, tau), 3), "\n", sep = "")
  cat("  propensity RMSE: ", round(sqrt(mean((propensity - p_carbo)^2)), 4),
      "\n", sep = "")
  cat("  propensity corr: ", round(cor(propensity, p_carbo), 3), "\n", sep = "")
}

# ---------------------------------------------------------------------------
# 8. Save
# ---------------------------------------------------------------------------

if (save_data) {
  ovarian       <- ovarian_new
  ovarian_truth <- ovarian_truth_new
  save(ovarian,       file = "data/ovarian.rda",       compress = "xz")
  save(ovarian_truth, file = "data/ovarian_truth.rda", compress = "xz")
  cat("\nWrote data/ovarian.rda and data/ovarian_truth.rda\n")
  cat("Update R/data-documentation.R, then run devtools::document().\n")
} else {
  cat("\nNothing written. Re-run with --save to overwrite data/*.rda\n")
}
if (!validate)
  cat("Re-run with --validate to refit the paper's models on this data.\n")
