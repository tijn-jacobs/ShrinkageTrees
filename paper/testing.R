remove.packages("ShrinkageTrees")
.rs.restartR()
setwd("/Users/tijnjacobs/Library/CloudStorage/OneDrive-VrijeUniversiteitAmsterdam/Documents/GitHub/ShrinkageTrees")
devtools::install()
devtools::load_all()
devtools::test()
devtools::document()
devtools::check()
?library(ShrinkageTrees)

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
                                    verbose = TRUE,
                                    n_chains = 3)

fit_double
summary(fit_double)                         




library(bayesplot)
library(ggplot2)

# ─────────────────────────────────────────────────────────────────────────────
# 1. ShrinkageTrees — continuous outcome, multiple chains
# ─────────────────────────────────────────────────────────────────────────────
set.seed(42)
n <- 150; p <- 8
X <- matrix(rnorm(n * p), n, p)
colnames(X) <- paste0("X", 1:p)
y <- X[, 1] * 2 + X[, 2] - 0.5 * X[, 3] + rnorm(n)

k  <- 0.5                         # BART k-factor
lh <- (max(y) - min(y)) / (k * 2 * sqrt(10))  # local_hp

fit_st <- ShrinkageTrees(
  y              = y,
  X_train        = X,
  prior_type     = "horseshoe",
  local_hp       = lh,
  global_hp      = lh,
  number_of_trees = 10,
  N_post         = 50,
  N_burn         = 20,
  store_posterior_sample = TRUE,
  n_chains       = 3,
  verbose        = TRUE
)

print(fit_st)
summary(fit_st)

# Plots
plot(fit_st, type = "trace")    # sigma traceplot — one line per chain
plot(fit_st, type = "density")  # overlaid posterior densities

# ─────────────────────────────────────────────────────────────────────────────
# 2. CausalShrinkageForest — continuous outcome, multiple chains
#    (single-chain first to verify the bug fix, then multi-chain)
# ─────────────────────────────────────────────────────────────────────────────
set.seed(123)
n   <- 200; p <- 6
X   <- matrix(rnorm(n * p), n, p)
colnames(X) <- paste0("X", 1:p)
w   <- rbinom(n, 1, 0.5)                          # treatment indicator
tau <- 1 + X[, 1]                                  # heterogeneous CATE
y   <- X[, 2] + tau * w + rnorm(n)

lh_c <- (max(y) - min(y)) / (0.5 * 2 * sqrt(10))

# -- single chain first (sanity check) --
fit_causal_1 <- CausalShrinkageForest(
  y                         = y,
  X_train_control           = X,
  X_train_treat             = X,
  treatment_indicator_train = w,
  prior_type_control        = "horseshoe",
  prior_type_treat          = "horseshoe",
  local_hp_control          = lh_c,
  global_hp_control         = lh_c,
  local_hp_treat            = lh_c,
  global_hp_treat           = lh_c,
  number_of_trees_control   = 10,
  number_of_trees_treat     = 10,
  N_post                    = 50,
  N_burn                    = 20,
  store_posterior_sample    = TRUE,
  n_chains                  = 1,
  verbose                   = TRUE
)
print(fit_causal_1)

# -- now multi-chain --
fit_causal <- CausalShrinkageForest(
  y                         = y,
  X_train_control           = X,
  X_train_treat             = X,
  treatment_indicator_train = w,
  prior_type_control        = "horseshoe",
  prior_type_treat          = "horseshoe",
  local_hp_control          = lh_c,
  global_hp_control         = lh_c,
  local_hp_treat            = lh_c,
  global_hp_treat           = lh_c,
  number_of_trees_control   = 10,
  number_of_trees_treat     = 10,
  N_post                    = 5000,
  N_burn                    = 20,
  store_posterior_sample    = TRUE,
  n_chains                  = 3,
  verbose                   = TRUE
)

print(fit_causal)
summary(fit_causal)

# Plots
plot(fit_causal, type = "trace")    # sigma mixing
plot(fit_causal, type = "density")  # sigma posterior
plot(fit_causal, type = "ate")      # ATE posterior density
plot(fit_causal, type = "cate") 
















library(ShrinkageTrees)

set.seed(42)
n <- 200; p <- 5
X <- matrix(rnorm(n * p), n, p)
W <- rbinom(n, 1, 0.5) 

# True CATE: heterogeneous, depends on X[,1] and X[,2]
tau_true <- 1.0 + 1.5 * X[, 1] - 0.8 * X[, 2]
mu_true  <- 2.0 * X[, 1] + X[, 3]
y <- mu_true + (W-0.5) * tau_true + rnorm(n, sd = 0.5)

# Propensity scores for adaptive coding
ps <- rep(0.5, n)  # known RCT design

# Common MCMC settings
args_common <- list(
  X_train_control = X, X_train_treat = X,
  treatment_indicator_train = W,
  number_of_trees = 50, N_post = 2000, N_burn = 2000,
  store_posterior_sample = TRUE, verbose = TRUE
)

# --- Centered ---
set.seed(1)
fit_cen <- do.call(CausalHorseForest, c(list(y = y, treatment_coding = "centered"), args_common))

# --- Binary ---
set.seed(1)
fit_bin <- do.call(CausalHorseForest, c(list(y = y, treatment_coding = "binary"), args_common))

# --- Adaptive ---
set.seed(1)
fit_ada <- do.call(CausalHorseForest, c(list(y = y, treatment_coding = "adaptive", propensity = ps), args_common))

# --- Invariant ---
set.seed(1)
fit_inv <- do.call(CausalHorseForest, c(list(y = y, treatment_coding = "invariant"), args_common))

# --- RMSE of CATE ---
rmse <- function(fit) sqrt(mean((fit$train_predictions_treat - tau_true)^2))

cat("\n=== CATE RMSE ===\n")
cat(sprintf("Centered:  %.4f\n", rmse(fit_cen)))
cat(sprintf("Binary:    %.4f\n", rmse(fit_bin)))
cat(sprintf("Adaptive:  %.4f\n", rmse(fit_ada)))
cat(sprintf("Invariant: %.4f\n", rmse(fit_inv)))

# --- ATE comparison ---
ate <- function(fit) mean(fit$train_predictions_treat)
cat("\n=== ATE (true:", round(mean(tau_true), 3), ") ===\n")
cat(sprintf("Centered:  %.4f\n", ate(fit_cen)))
cat(sprintf("Binary:    %.4f\n", ate(fit_bin)))
cat(sprintf("Adaptive:  %.4f\n", ate(fit_ada)))
cat(sprintf("Invariant: %.4f\n", ate(fit_inv)))

# --- Invariant: inspect b0/b1 ---
cat("\n=== Invariant coding parameters ===\n")
cat(sprintf("b0: mean = %.3f, sd = %.3f\n", mean(fit_inv$b0), sd(fit_inv$b0)))
cat(sprintf("b1: mean = %.3f, sd = %.3f\n", mean(fit_inv$b1), sd(fit_inv$b1)))
cat(sprintf("b1 - b0: mean = %.3f\n", mean(fit_inv$b1 - fit_inv$b0)))