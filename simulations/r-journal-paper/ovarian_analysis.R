###############################################################################
# Worked example from the ShrinkageTrees R Journal paper, Sections 2 and 4.
#
# The analysis chunks of revisions/revision_v1/ShrinkageTrees.Rmd, extracted
# verbatim and in the order the manuscript runs them, so the numbers here are
# the numbers reported there.
#
# One deviation: the two tables are rendered with knitr::kable() in the
# manuscript and printed as plain data frames here. The quantities are the
# same; only the formatting differs.
#
# Requires ShrinkageTrees >= 2.1.0.
###############################################################################

library(ShrinkageTrees)
library(survival)
library(ggplot2)

stopifnot(utils::packageVersion("ShrinkageTrees") >= "2.1.0")

# ---- Data and design matrix ----------------------------------------------------

library(ShrinkageTrees)
data("ovarian")
data("ovarian_truth")
time   <- ovarian$OS_time / 30.44  # days to months
status <- ovarian$OS_event
X      <- as.matrix(ovarian[, setdiff(names(ovarian),
                                      c("OS_time", "OS_event"))])
fit    <- HorseTrees(y = time, status = status, X_train = X,
                     outcome_type = "right-censored", number_of_trees = 200,
                     verbose = FALSE)

treatment <- ovarian$treatment

# ---- Train / test split --------------------------------------------------------

set.seed(42)
test_idx  <- sample(seq_len(nrow(X)), size = round(0.3 * nrow(X)))
train_idx <- setdiff(seq_len(nrow(X)), test_idx)

# ---- Section 4.1  Standard BART baseline ---------------------------------------

fit_bart <- SurvivalBART(
  time            = time[train_idx],
  status          = status[train_idx],
  X_train         = X[train_idx, ],
  X_test          = X[test_idx, ],
  timescale       = "time",
  number_of_trees = 200,
  N_post          = 5000,
  N_burn          = 5000,
  n_chains        = 4,
  verbose         = FALSE
)

print(fit_bart)

library(survival)
c_bart_in <- concordance(Surv(time[train_idx], status[train_idx]) ~
                           fit_bart$train_predictions)$concordance

plot(fit_bart, type = "survival", km = TRUE)

# ---- Section 4.1  Horseshoe Forest ---------------------------------------------

fit_horse <- HorseTrees(
  y               = time[train_idx],
  status          = status[train_idx],
  X_train         = X[train_idx, ],
  X_test          = X[test_idx, ],
  outcome_type    = "right-censored",
  timescale       = "time",
  number_of_trees = 200,
  k               = 1,
  N_post          = 5000,
  N_burn          = 5000,
  n_chains        = 4,
  verbose         = FALSE
)

# ---- Section 4.1  Prediction metrics -------------------------------------------

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

f_true  <- ovarian_truth$f
train   <- list(`SurvivalBART()` = log(fit_bart$train_predictions),
                `HorseTrees()`   = log(fit_horse$train_predictions))
test    <- list(`SurvivalBART()` = log(fit_bart$test_predictions),
                `HorseTrees()`   = log(fit_horse$test_predictions))
train_s <- list(`SurvivalBART()` = log(fit_bart$train_predictions_sample),
                `HorseTrees()`   = log(fit_horse$train_predictions_sample))
test_s  <- list(`SurvivalBART()` = log(fit_bart$test_predictions_sample),
                `HorseTrees()`   = log(fit_horse$test_predictions_sample))

mt <- data.frame(
  Model = names(train),
  ci_in  = vapply(train, cidx, 0, i = train_idx),
  ci_out = vapply(test,  cidx, 0, i = test_idx),
  rm_in  = vapply(train, rmse, 0, b = f_true[train_idx]),
  rm_out = vapply(test,  rmse, 0, b = f_true[test_idx]),
  cv_in  = vapply(train_s, covg, 0, b = f_true[train_idx]),
  cv_out = vapply(test_s,  covg, 0, b = f_true[test_idx]),
  wd_in  = vapply(train_s, widt, 0),
  wd_out = vapply(test_s,  widt, 0),
  row.names = NULL
)

mt

plot(fit_horse, type = "trace")
plot(fit_horse, type = "density")

# ---- Section 4.1  Posterior survival curve -------------------------------------

idx  <- sample(test_idx, 1)
pred <- predict(fit_horse, newdata = X[idx, , drop = FALSE])

plot(pred, type = "survival", obs = 1)

# ---- Section 4.2  Propensity scores --------------------------------------------

X_cov <- X[, setdiff(colnames(X), "treatment")]
ps_fit <- HorseTrees(
  y               = treatment,
  X_train         = X_cov,
  outcome_type    = "binary",
  number_of_trees = 200,
  k               = 1,
  N_post          = 5000,
  N_burn          = 5000,
  verbose         = FALSE
)
propensity <- ps_fit$train_predictions
X_control  <- cbind(propensity = propensity, X_cov)

# ---- Section 4.2  Causal forest ------------------------------------------------

fit_causal <- CausalHorseForest(
  y                         = log(time),
  status                    = status,
  X_train_control           = X_control,
  X_train_treat             = X_cov,
  treatment_indicator_train = treatment,
  outcome_type              = "right-censored",
  timescale                 = "log",
  number_of_trees           = 200,
  k                         = 1.5,
  N_post                    = 5000,
  N_burn                    = 5000,
  store_posterior_sample    = TRUE,
  n_chains                  = 4,
  verbose                   = FALSE
)
causal_summary <- summary(fit_causal)
causal_summary

plot(fit_causal, type = "ate")
plot(fit_causal, type = "cate")

# ---- Section 4.2  Treatment effect against the truth ---------------------------

te   <- causal_summary$treatment_effect
cate <- fit_causal$train_predictions_treat
ci   <- apply(fit_causal$train_predictions_sample_treat, 2, quantile,
              c(0.025, 0.975))
tau_true <- ovarian_truth$tau

ate_true   <- mean(tau_true)
ate_cov    <- ate_true >= te$ate_lower && ate_true <= te$ate_upper
cate_rmse  <- sqrt(mean((cate - tau_true)^2))
cate_cov   <- mean(tau_true >= ci[1, ] & tau_true <= ci[2, ])
cate_width <- mean(ci[2, ] - ci[1, ])

data.frame(`CATE RMSE`  = round(cate_rmse, 3),
           `CATE cov`   = round(cate_cov, 3),
           `CATE width` = round(cate_width, 3),
           check.names = FALSE, row.names = NULL)

# ---- Section 4.3  Posterior projection -----------------------------------------

# Project the CATE surface onto an additive model with selected covariates
proj_tau <- posterior_projection(fit_causal, target = "treatment",
                                 family = "additive",
                                 covariates = c("figo_stage", "age"))
proj_tau                            # coefficient table, summary R^2
plot(proj_tau, type = "effects")    # per-covariate curves, credible bands

# Sparse summary of the posterior mean CATE surface
posterior_projection(fit_causal, target = "treatment",
                     family = "linear", penalty = "lasso")

