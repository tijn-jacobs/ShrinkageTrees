# ============================================================
# Continuous outcome simulation: uniform evaluators + metrics
# Methods: BART, HorseTrees (HS), DART, SoftBART, RF, Elastic Net, MARS
# Notes:
# - Uses dt$X and dt$y consistently (your current code mixes X vs X_train).
# - For Bayesian tree methods, returns BOTH point + posterior metrics.
# - For non-Bayesian methods, returns point metrics and NA for posterior metrics.
# ============================================================

suppressPackageStartupMessages({
  library(ShrinkageTrees)
  library(foreach)
  library(doParallel)
})

# -------------------------------
# DGP
# -------------------------------
data_gen_continuous <- function(n, p, sigma, s = 0.1,
                                nonlinear = TRUE, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  X <- matrix(runif(n * p), n, p)
  
  beta <- rnorm(p, 0, 1) * rbinom(p, 1, s)
  mu_linear <- as.vector(X %*% beta)
  
  if (nonlinear) {
    f <- 10 * sin(pi * X[, 1] * X[, 2]) +
      20 * (X[, 3] - 0.5)^2 +
      10 * X[, 4] +
      5  * X[, 5]
  } else {
    f <- mu_linear
  }
  
  y <- f + rnorm(n, 0, sigma)
  
  list(X = X, y = y, f_true = f)
}

# -------------------------------
# Metrics (point + posterior)
# -------------------------------
compute_point_metrics <- function(y_true, y_pred) {
  rmse <- sqrt(mean((y_true - y_pred)^2))
  mae  <- mean(abs(y_true - y_pred))
  r2   <- 1 - sum((y_true - y_pred)^2) / sum((y_true - mean(y_true))^2)
  
  list(RMSE = rmse, MAE = mae, R2 = r2)
}

compute_posterior_metrics <- function(y_true, y_pred_mean, y_pred_samps,
                                      alpha = 0.05) {
  ci <- apply(y_pred_samps, 2, quantile, probs = c(alpha / 2, 1 - alpha / 2))
  coverage  <- mean(ci[1, ] <= y_true & ci[2, ] >= y_true)
  ci_length <- mean(ci[2, ] - ci[1, ])
  
  # CRPS (empirical, using posterior predictive draws)
  # E|X - y| - 0.5 E|X - X'|
  crps <- mean(
    colMeans(abs(y_pred_samps - rep(y_true, each = nrow(y_pred_samps)))) -
      0.5 * apply(y_pred_samps, 2, function(z) {
        mean(abs(outer(z, z, "-")))
      })
  )
  
  c(
    compute_point_metrics(y_true, y_pred_mean),
    list(Coverage = coverage, CI_Length = ci_length, CRPS = crps)
  )
}

as_row <- function(metrics, method) {
  df <- as.data.frame(metrics, optional = TRUE)
  df$Method <- method
  df
}

# -------------------------------
# Helpers
# -------------------------------
default_local_hp <- function(y, m = 200) (max(y) - min(y)) / (2 * sqrt(m))

fill_posterior_nas <- function(point_metrics) {
  c(point_metrics, list(Coverage = NA_real_, CI_Length = NA_real_, CRPS = NA_real_))
}

# ============================================================
# Evaluators
# ============================================================

evaluate_BART <- function(dt, m = 200, N_post = 1000, N_burn = 1000, ...) {
  fit <- ShrinkageTrees::ShrinkageTrees(
    y = dt$y,
    X_train = dt$X,
    outcome_type = "continuous",
    prior_type = "standard",
    number_of_trees = m,
    local_hp = 1/sqrt(m),#default_local_hp(dt$y, m),
    store_posterior_sample = TRUE,
    N_post = N_post,
    N_burn = N_burn,
    verbose = FALSE,
    ...
  )
  
  compute_posterior_metrics(
    y_true = dt$f_true,
    y_pred_mean = fit$train_predictions,
    y_pred_samps = fit$train_predictions_sample
  )
}

evaluate_HorseTrees <- function(dt, m = 200,
                                k_local = 0.1, k_global = 0.1,
                                N_post = 1000, N_burn = 1000, ...) {
  fit <- ShrinkageTrees::ShrinkageTrees(
    y = dt$y,
    X_train = dt$X,
    outcome_type = "continuous",
    prior_type = "horseshoe",
    number_of_trees = m,
    local_hp = k_local / sqrt(m),
    global_hp = k_global / sqrt(m),
    store_posterior_sample = TRUE,
    N_post = N_post,
    N_burn = N_burn,
    verbose = FALSE,
    ...
  )
  
  compute_posterior_metrics(
    y_true = dt$f_true,
    y_pred_mean = fit$train_predictions,
    y_pred_samps = fit$train_predictions_sample
  )
}

evaluate_DART <- function(dt, m = 200, N_post = 1000, N_burn = 1000, ...) {
  # Assumption: your corrected package exposes a DART-like prior_type.
  # If your package uses a different flag/name, change prior_type below.
  fit <- ShrinkageTrees::ShrinkageTrees(
    y = dt$y,
    X_train = dt$X,
    outcome_type = "continuous",
    prior_type = "dart",
    number_of_trees = m,
    local_hp = default_local_hp(dt$y, m),
    store_posterior_sample = TRUE,
    N_post = N_post,
    N_burn = N_burn,
    verbose = FALSE,
    ...
  )
  
  compute_posterior_metrics(
    y_true = dt$f_true,
    y_pred_mean = fit$train_predictions,
    y_pred_samps = fit$train_predictions_sample
  )
}

evaluate_SoftBART <- function(dt, m = 200, N_post = 1000, N_burn = 1000, ...) {
  # If SoftBART is a different entry point in your package, swap it here.
  # This keeps the interface uniform.
  fit <- ShrinkageTrees::ShrinkageTrees(
    y = dt$y,
    X_train = dt$X,
    outcome_type = "continuous",
    prior_type = "softbart",
    number_of_trees = m,
    local_hp = default_local_hp(dt$y, m),
    store_posterior_sample = TRUE,
    N_post = N_post,
    N_burn = N_burn,
    verbose = FALSE,
    ...
  )
  
  compute_posterior_metrics(
    y_true = dt$f_true,
    y_pred_mean = fit$train_predictions,
    y_pred_samps = fit$train_predictions_sample
  )
}

evaluate_RF <- function(dt, ntree = 500, ...) {
  if (!requireNamespace("randomForest", quietly = TRUE)) {
    stop("Package 'randomForest' required.")
  }
  fit <- randomForest::randomForest(x = dt$X, y = dt$y, ntree = ntree, ...)
  y_pred <- as.numeric(predict(fit, dt$X))
  fill_posterior_nas(compute_point_metrics(dt$f_true, y_pred))
}

evaluate_ElasticNet <- function(dt, alpha = 0.5, ...) {
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("Package 'glmnet' required.")
  }
  fit <- glmnet::cv.glmnet(
    x = dt$X, y = dt$y, alpha = alpha,
    standardize = TRUE, ...
  )
  y_pred <- as.numeric(predict(fit, newx = dt$X, s = "lambda.min"))
  fill_posterior_nas(compute_point_metrics(dt$f_true, y_pred))
}

evaluate_MARS <- function(dt, ...) {
  if (!requireNamespace("earth", quietly = TRUE)) {
    stop("Package 'earth' required (MARS).")
  }
  fit <- earth::earth(x = dt$X, y = dt$y, ...)
  y_pred <- as.numeric(predict(fit, dt$X))
  fill_posterior_nas(compute_point_metrics(dt$f_true, y_pred))
}

# ============================================================
# Single simulation
# ============================================================
run_single_simulation_regression <- function(n, p, sigma, s,
                                             nonlinear,
                                             N_post = 1000, N_burn = 1000,
                                             seed = NULL) {
  dt <- data_gen_continuous(
    n = n, p = p, sigma = sigma, s = s,
    nonlinear = nonlinear, seed = seed
  )
  
  res <- list(
    BART      = evaluate_BART(dt, N_post = N_post, N_burn = N_burn),
    HorseTrees = evaluate_HorseTrees(dt, N_post = N_post, N_burn = N_burn),
    DART      = evaluate_DART(dt, N_post = N_post, N_burn = N_burn),
    SoftBART  = evaluate_SoftBART(dt, N_post = N_post, N_burn = N_burn),
    RF        = evaluate_RF(dt),
    ElasticNet = evaluate_ElasticNet(dt, alpha = 0.5),
    MARS      = evaluate_MARS(dt)
  )
  
  out <- do.call(rbind, Map(as_row, res, names(res)))
  rownames(out) <- NULL
  out$n <- n
  out$p <- p
  out$sigma <- sigma
  out$sparsity <- s
  out$nonlinear <- nonlinear
  out$seed <- if (is.null(seed)) NA_integer_ else seed
  out
}

# ============================================================
# Parallel study
# ============================================================
run_simulation_study <- function(M, n, p, sigma, s, nonlinear,
                                 N_post = 1000, N_burn = 1000,
                                 n_cores = 2) {
  doParallel::registerDoParallel(cores = n_cores)
  
  foreach(i = 1:M, .combine = rbind,
          .packages = c("ShrinkageTrees", "randomForest", "glmnet", "earth")) %dopar% {
            run_single_simulation_regression(
              n = n, p = p, sigma = sigma, s = s,
              nonlinear = nonlinear,
              N_post = N_post, N_burn = N_burn,
              seed = i
            )
          }
}

# Example:
# res <- run_simulation_study(M=50, n=200, p=1000, sigma=1, s=0.1,
#                             nonlinear=TRUE, N_post=2000, N_burn=1000,
#                             n_cores=8)
# aggregate(. ~ Method, res[,c("Method","RMSE","MAE","R2","Coverage","CI_Length","CRPS")], mean, na.rm=TRUE)












library(ShrinkageTrees)
library(randomForest)
library(glmnet)

set.seed(1)

# --------------------------------------------------
# Generate one dataset
# --------------------------------------------------
dt <- data_gen_continuous(
  n = 200,
  p = 50,
  sigma = 1,
  s = 0.1,
  nonlinear = TRUE,
  seed = 1
)

# For consistency with your evaluators
dt$X_train <- dt$X

# --------------------------------------------------
# Run each method separately
# --------------------------------------------------

res_BART <- evaluate_BART(
  dt,
  N_post = 500,
  N_burn = 500
)

res_HorseTrees <- evaluate_HorseTrees(
  dt = dt,
  local_hp = 0.1,
  global_hp = 0.1,
  N_post = 500,
  N_burn = 500
)

res_DART <- evaluate_DART(
  dt = dt,
  N_post = 500,
  N_burn = 500
)

res_RF <- evaluate_RF(
  dt = dt,
  ntree = 500
)

res_EN <- evaluate_ElasticNet(
  dt = dt,
  alpha = 0.5
)

# --------------------------------------------------
# Collect results in a clean table
# --------------------------------------------------
results <- rbind(
  BART        = unlist(res_BART),
  HorseTrees = unlist(res_HorseTrees),
  DART        = unlist(res_DART),
  RF          = unlist(res_RF),
  ElasticNet  = unlist(res_EN)
)

print(round(results, 3))



