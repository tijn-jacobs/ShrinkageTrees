remove.packages("ShrinkageTrees")
.rs.restartR()
setwd("/Users/tijnjacobs/Library/CloudStorage/OneDrive-VrijeUniversiteitAmsterdam/Documents/GitHub/ShrinkageTrees")
devtools::install()
devtools::load_all()
library(ShrinkageTrees)

n <- 50
p <- 3
X <- matrix(runif(n * p), ncol = p)
X_test <- matrix(runif(n * p), ncol = p)
y <- X[, 1] + rnorm(n)

# Fit ShrinkageTrees with standard horseshoe prior
fit_horseshoe <- ShrinkageTrees(y = y,
                                X_train = X,
                                X_test = NULL,
                                outcome_type = "continuous",
                                number_of_trees = 5,
                                prior_type = "horseshoe",
                                local_hp = 0.1 / sqrt(5),
                                global_hp = 0.1 / sqrt(5),
                                N_post = 1,
                                N_burn = 5,
                                store_posterior_sample = TRUE,
                                verbose = FALSE,
                                n_chains = 5)

fit_horseshoe
summary(fit_horseshoe)

# Predict on new data
pred <- predict(fit_horseshoe, newdata = X_test)
str(pred)










# Example: Continuous outcome, homogenuous treatment effect, two priors
n <- 50
p <- 3
X <- matrix(runif(n * p), ncol = p)
X_treat <- X_control <- X
treat <- rbinom(n, 1, X[,1])
tau <- 2
y <- X[, 1] + (0.5 - treat) * tau + rnorm(n)

# Fit a standard Causal Horseshoe Forest
fit_double <- CausalShrinkageForest(y = y,
                                    X_train_control = X_control,
                                    X_train_treat = X_treat,
                                    treatment_indicator_train = treat,
                                    outcome_type = "continuous",
                                    number_of_trees_treat = 5,
                                    number_of_trees_control = 5,
                                    prior_type_control = "horseshoe",
                                    prior_type_treat = "horseshoe",
                                    local_hp_control = 0.1/sqrt(5),
                                    local_hp_treat = 0.1/sqrt(5),
                                    global_hp_control = 0.1/sqrt(5),
                                    global_hp_treat = 0.1/sqrt(5),
                                    N_post = 100,
                                    N_burn = 5,
                                    store_posterior_sample = TRUE,
                                    verbose = FALSE)

fit_double
summary(fit_double)                         





















# ============================================================
# Test script for SurvivalBART, SurvivalDART, SurvivalBCF,
# SurvivalShrinkageBCF
# ============================================================

library(survival)

# Shared data
set.seed(42)
n <- 50
p <- 3
X <- matrix(runif(n * p), ncol = p)
X_test <- matrix(runif(n * p), ncol = p)

# Simulate right-censored survival times
true_time <- exp(X[, 1] + rnorm(n, sd = 0.5))
cens_time <- runif(n, 0.5, 3)
time   <- pmin(true_time, cens_time)
status <- as.integer(true_time <= cens_time)
treat  <- rbinom(n, 1, X[, 1])

# ============================================================
# 1. SurvivalBART
# ============================================================
fit_sbart <- SurvivalBART(
  time             = time,
  status           = status,
  X_train          = X,
  X_test           = X_test,
  timescale        = "time",
  number_of_trees  = 5,
  k                = 2.0,
  N_post           = 10,
  N_burn           = 5,
  verbose          = FALSE
)
fit_sbart
summary(fit_sbart)

# ============================================================
# 2. SurvivalDART
# ============================================================
fit_sdart <- SurvivalDART(
  time             = time,
  status           = status,
  X_train          = X,
  X_test           = X_test,
  timescale        = "time",
  number_of_trees  = 5,
  a_dirichlet      = 0.5,
  b_dirichlet      = 1.0,
  k                = 2.0,
  N_post           = 10,
  N_burn           = 5,
  verbose          = FALSE
)
fit_sdart
summary(fit_sdart)

# ============================================================
# 3. SurvivalBCF
# ============================================================
fit_sbcf <- SurvivalBCF(
  time                    = time,
  status                  = status,
  X_train                 = X,
  treatment               = treat,
  timescale               = "time",
  propensity              = NULL,
  number_of_trees_control = 5,
  number_of_trees_treat   = 5,
  N_post                  = 10,
  N_burn                  = 5,
  verbose                 = FALSE
)
fit_sbcf
summary(fit_sbcf)

# ============================================================
# 4. SurvivalShrinkageBCF
# ============================================================
fit_ssbcf <- SurvivalShrinkageBCF(
  time                    = time,
  status                  = status,
  X_train                 = X,
  treatment               = treat,
  timescale               = "time",
  propensity              = NULL,
  a_dir                   = 0.5,
  b_dir                   = 1.0,
  number_of_trees_control = 5,
  number_of_trees_treat   = 5,
  N_post                  = 10,
  N_burn                  = 5,
  verbose                 = FALSE
)
fit_ssbcf
summary(fit_ssbcf)              
