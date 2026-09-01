###############################################################################
# Worked example from the ShrinkageTrees R Journal paper, Sections 2 and 4.
# The analysis chunks in the order the manuscript runs them, at REDUCED MCMC
# settings for runtime (N_post = N_burn = 1000 here; the manuscript chunks use
# 5000/5000). Numbers here therefore approximate the published ones.
# Requires ShrinkageTrees >= 2.1.0.
###############################################################################

library(ShrinkageTrees)
library(survival)
library(ggplot2)

# ---- Data and design matrix -------------------------------------------------

data("ovarian")
data("ovarian_truth")

time      <- ovarian$OS_time / 30.44  # days to months
status    <- ovarian$OS_event
treatment <- ovarian$treatment

# Prediction models treat treatment as an ordinary covariate; the causal models
# take it separately, so it is excluded from their design.
X_pred <- as.matrix(ovarian[, setdiff(names(ovarian),
                                      c("OS_time", "OS_event"))])
X      <- as.matrix(ovarian[, setdiff(names(ovarian),
                                      c("OS_time", "OS_event", "treatment"))])

fit <- HorseTrees(y = time, status = status, X_train = X_pred,
                  outcome_type = "right-censored", number_of_trees = 200,
                  verbose = FALSE)

# ---- Section 4.1  Standard BART baseline ------------------------------------

set.seed(42)
test_idx  <- sample(seq_len(nrow(X_pred)), size = round(0.3 * nrow(X_pred)))
train_idx <- setdiff(seq_len(nrow(X_pred)), test_idx)

fit_bart <- SurvivalBART(
  time            = time[train_idx],
  status          = status[train_idx],
  X_train         = X_pred[train_idx, ],
  X_test          = X_pred[test_idx, ],
  timescale       = "time",
  number_of_trees = 200,
  N_post          = 1000,
  N_burn          = 1000,
  n_chains        = 4,
  verbose         = FALSE
)

print(fit_bart)
plot(fit_bart, type = "survival", km = TRUE)

# ---- Section 4.1  Horseshoe Forest ------------------------------------------

fit_horse <- HorseTrees(
  y               = time[train_idx],
  status          = status[train_idx],
  X_train         = X_pred[train_idx, ],
  X_test          = X_pred[test_idx, ],
  outcome_type    = "right-censored",
  timescale       = "time",
  number_of_trees = 200,
  k               = 1,
  N_post          = 1000,
  N_burn          = 1000,
  n_chains        = 4,
  verbose         = FALSE
)

plot(fit_horse, type = "trace")
plot(fit_horse, type = "density")

# ---- Section 4.1  Prediction metrics ----------------------------------------

cidx <- function(p, i) concordance(Surv(time[i], status[i]) ~ p)$concordance
rmse <- function(a, b) sqrt(mean((a - b)^2))
covg <- function(s, b) {
  ci <- apply(s, 2, quantile, c(0.025, 0.975))
  mean(b >= ci[1, ] & b <= ci[2, ])
}
widt <- function(s) {
  ci <- apply(s, 2, quantile, c(0.025, 0.975))
  mean(ci[2, ] - ci[1, ])
}

# Coverage is the share of patients whose true value lies in their own 95%
# interval. ovarian_truth$f is the regression function the models estimate.
f_true  <- ovarian_truth$f
train   <- list(SurvivalBART = log(fit_bart$train_predictions),
                HorseTrees   = log(fit_horse$train_predictions))
test    <- list(SurvivalBART = log(fit_bart$test_predictions),
                HorseTrees   = log(fit_horse$test_predictions))
train_s <- list(SurvivalBART = log(fit_bart$train_predictions_sample),
                HorseTrees   = log(fit_horse$train_predictions_sample))
test_s  <- list(SurvivalBART = log(fit_bart$test_predictions_sample),
                HorseTrees   = log(fit_horse$test_predictions_sample))

data.frame(
  Model      = names(train),
  `C in`     = round(vapply(train, cidx, 0, i = train_idx), 3),
  `C out`    = round(vapply(test,  cidx, 0, i = test_idx),  3),
  `RMSE in`  = round(vapply(train, rmse, 0, b = f_true[train_idx]), 3),
  `RMSE out` = round(vapply(test,  rmse, 0, b = f_true[test_idx]),  3),
  `Cov in`    = round(vapply(train_s, covg, 0, b = f_true[train_idx]), 3),
  `Cov out`   = round(vapply(test_s,  covg, 0, b = f_true[test_idx]),  3),
  `Width in`  = round(vapply(train_s, widt, 0), 3),
  `Width out` = round(vapply(test_s,  widt, 0), 3),
  check.names = FALSE, row.names = NULL
)

# ---- Section 4.2  Propensity scores -----------------------------------------

ps_fit <- HorseTrees(
  y               = treatment,
  X_train         = X,
  outcome_type    = "binary",
  number_of_trees = 200,
  k               = 1,
  N_post          = 1000,
  N_burn          = 1000,
  verbose         = FALSE
)
X_control <- cbind(propensity = ps_fit$train_predictions, X)

# ---- Section 4.2  Causal forest ---------------------------------------------

fit_causal <- CausalHorseForest(
  y                         = log(time),
  status                    = status,
  X_train_control           = X_control,
  X_train_treat             = X,
  treatment_indicator_train = treatment,
  outcome_type              = "right-censored",
  timescale                 = "log",
  number_of_trees           = 200,
  k                         = 1.5,
  N_post                    = 1000,
  N_burn                    = 1000,
  store_posterior_sample    = TRUE,
  n_chains                  = 4,
  verbose                   = FALSE
)

causal_summary <- summary(fit_causal)
causal_summary
plot(fit_causal, type = "ate")
plot(fit_causal, type = "cate")

# Individual interval coverage: the share of patients whose true CATE falls in
# their own 95% credible interval. Nominal is 0.95.
cate <- fit_causal$train_predictions_treat
ci   <- apply(fit_causal$train_predictions_sample_treat, 2, quantile,
              c(0.025, 0.975))

# summary() reports the population ATE by default, via Bayesian bootstrap.
te       <- causal_summary$treatment_effect
tau_true <- ovarian_truth$tau

# PATE against the truth, reported in the text rather than the table.
c(pate = round(te$ate, 3), true = round(mean(tau_true), 3),
  lower = round(te$ate_lower, 3), upper = round(te$ate_upper, 3))

data.frame(
  Model        = "CausalHorseForest",
  `CATE RMSE`  = round(rmse(cate, tau_true), 3),
  `CATE cov`   = round(mean(tau_true >= ci[1, ] & tau_true <= ci[2, ]), 3),
  `CATE width` = round(mean(ci[2, ] - ci[1, ]), 3),
  check.names = FALSE, row.names = NULL
)

# ---- Section 4.3  Posterior projection --------------------------------------

# tau depends on figo_stage, age, and the 3rd and 7th gene, and on nothing
# else. Projecting onto exactly that set should recover it.
genes  <- grep("^ENSG", names(ovarian), value = TRUE)
active <- c("figo_stage", "age", genes[c(3, 7)])

proj <- posterior_projection(fit_causal, target = "treatment",
                             family = "additive", covariates = active)
proj
plot(proj, type = "effects")
