test_that("CausalHorseForest works for continuous outcome", {
  
  # Generate data
  n <- 50
  p <- 3
  X_control <- matrix(runif(n * p), ncol = p)
  X_treat   <- matrix(runif(n * p), ncol = p)
  treatment <- rbinom(n, 1, 0.5)
  
  # True treatment effect
  tau <- 2
  
  # Outcome: baseline + treatment effect + noise
  y <- X_control[, 1] + treatment * tau + rnorm(n)
  
  # Test data
  X_test_control <- matrix(runif(n * p), ncol = p)
  X_test_treat   <- matrix(runif(n * p), ncol = p)
  treatment_test <- rbinom(n, 1, 0.5)
  
  # Fit the model
  fit <- CausalHorseForest(
    y = y,
    X_train_control = X_control,
    X_train_treat = X_treat,
    treatment_indicator_train = treatment,
    X_test_control = X_test_control,
    X_test_treat = X_test_treat,
    treatment_indicator_test = treatment_test,
    outcome_type = "continuous",
    number_of_trees = 5,
    N_post = 10,
    N_burn = 5,
    store_posterior_sample = TRUE,
    verbose = FALSE,
    seed = 1
  )
  
  # --- Basic list checks ---
  expect_type(fit, "list")
  
  # Check scalar outputs
  expect_length(fit$acceptance_ratio_control, 1)
  expect_length(fit$acceptance_ratio_treat, 1)
  expect_true(all(fit$acceptance_ratio_control >= 0 & fit$acceptance_ratio_control <= 1))
  expect_true(all(fit$acceptance_ratio_treat >= 0 & fit$acceptance_ratio_treat <= 1))
  
  # --- Check main predictions ---
  expect_length(fit$train_predictions, n)
  expect_length(fit$test_predictions, n)
  expect_length(fit$train_predictions_control, n)
  expect_length(fit$test_predictions_control, n)
  expect_length(fit$train_predictions_treat, n)
  expect_length(fit$test_predictions_treat, n)
  
  # --- Check sigma vector ---
  expect_length(fit$sigma, 10 + 5)
  expect_true(all(fit$sigma > 0))
  
  # --- Posterior samples (control) ---
  expect_equal(dim(fit$train_predictions_sample_control), c(10, n))
  expect_equal(dim(fit$test_predictions_sample_control), c(10, n))
  
  # --- Posterior samples (treat) ---
  expect_equal(dim(fit$train_predictions_sample_treat), c(10, n))
  expect_equal(dim(fit$test_predictions_sample_treat), c(10, n))
  
  # --- Numerical sanity checks ---
  expect_false(any(is.na(unlist(fit))))
  expect_false(any(is.nan(unlist(fit))))
  expect_true(all(is.finite(unlist(fit))))
  
  # Predictions are not all zero
  expect_true(sd(fit$train_predictions) > 0)
  expect_true(sd(fit$test_predictions) > 0)
  expect_true(sd(fit$train_predictions_control) > 0)
  expect_true(sd(fit$test_predictions_control) > 0)
  expect_true(sd(fit$train_predictions_treat) > 0)
  expect_true(sd(fit$test_predictions_treat) > 0)
  
  # Posterior samples non-constant
  expect_true(sd(fit$train_predictions_sample_control) > 0)
  expect_true(sd(fit$test_predictions_sample_control) > 0)
  expect_true(sd(fit$train_predictions_sample_treat) > 0)
  expect_true(sd(fit$test_predictions_sample_treat) > 0)
  
  # --- Check reproducibility ---
  fit2 <- CausalHorseForest(
    y = y,
    X_train_control = X_control,
    X_train_treat = X_treat,
    treatment_indicator_train = treatment,
    outcome_type = "continuous",
    number_of_trees = 5,
    N_post = 10,
    N_burn = 5,
    store_posterior_sample = TRUE,
    verbose = FALSE,
    seed = 1
  )
  expect_equal(fit$train_predictions, fit2$train_predictions)
})


test_that("CausalHorseForest works for survival outcome", {
  
  # Generate data
  n <- 50
  p <- 3
  X_control <- matrix(runif(n * p), ncol = p)
  X_treat   <- matrix(runif(n * p), ncol = p)
  treatment <- rbinom(n, 1, 0.5)
  
  # True treatment effect on log-time (AFT interpretation)
  tau <- -0.5  # Negative tau => accelerates survival (protective)
  
  # Baseline log survival time
  log_T0 <- X_control[, 1] + rnorm(n)
  
  # Treatment effect
  log_T <- log_T0 + treatment * tau
  
  # Transform to actual survival time
  Ti <- exp(log_T)
  
  # Administrative censoring times (e.g., follow-up cutoff)
  C <- rexp(n, rate = 0.1)  # example: random censoring times
  
  # Observed time and status indicator
  time <- pmin(Ti, C)
  status <- as.integer(T <= C)  # 1 = event (uncensored), 0 = censored
  
  # Test data
  X_test_control <- matrix(runif(n * p), ncol = p)
  X_test_treat   <- matrix(runif(n * p), ncol = p)
  treatment_test <- rbinom(n, 1, 0.5)
  
  # Fit the model
  fit <- CausalHorseForest(
    y = time,
    status = status,
    X_train_control = X_control,
    X_train_treat = X_treat,
    treatment_indicator_train = treatment,
    X_test_control = X_test_control,
    X_test_treat = X_test_treat,
    treatment_indicator_test = treatment_test,
    outcome_type = "right-censored",
    number_of_trees = 5,
    N_post = 10,
    N_burn = 5,
    store_posterior_sample = TRUE,
    verbose = FALSE,
    seed = 1
  )
  
  # --- Basic list checks ---
  expect_type(fit, "list")
  
  # Check scalar outputs
  expect_length(fit$acceptance_ratio_control, 1)
  expect_length(fit$acceptance_ratio_treat, 1)
  expect_true(all(fit$acceptance_ratio_control >= 0 & fit$acceptance_ratio_control <= 1))
  expect_true(all(fit$acceptance_ratio_treat >= 0 & fit$acceptance_ratio_treat <= 1))
  
  # --- Check main predictions ---
  expect_length(fit$train_predictions, n)
  expect_length(fit$test_predictions, n)
  expect_length(fit$train_predictions_control, n)
  expect_length(fit$test_predictions_control, n)
  expect_length(fit$train_predictions_treat, n)
  expect_length(fit$test_predictions_treat, n)
  
  # --- Check sigma vector ---
  expect_length(fit$sigma, 10 + 5)
  expect_true(all(fit$sigma > 0))
  
  # --- Posterior samples (control) ---
  expect_equal(dim(fit$train_predictions_sample_control), c(10, n))
  expect_equal(dim(fit$test_predictions_sample_control), c(10, n))
  
  # --- Posterior samples (treat) ---
  expect_equal(dim(fit$train_predictions_sample_treat), c(10, n))
  expect_equal(dim(fit$test_predictions_sample_treat), c(10, n))
  
  # --- Numerical sanity checks ---
  expect_false(any(is.na(unlist(fit))))
  expect_false(any(is.nan(unlist(fit))))
  expect_true(all(is.finite(unlist(fit))))
  
  # Predictions are not all zero
  expect_true(sd(fit$train_predictions) > 0)
  expect_true(sd(fit$test_predictions) > 0)
  expect_true(sd(fit$train_predictions_control) > 0)
  expect_true(sd(fit$test_predictions_control) > 0)
  expect_true(sd(fit$train_predictions_treat) > 0)
  expect_true(sd(fit$test_predictions_treat) > 0)
  
  # Posterior samples non-constant
  expect_true(sd(fit$train_predictions_sample_control) > 0)
  expect_true(sd(fit$test_predictions_sample_control) > 0)
  expect_true(sd(fit$train_predictions_sample_treat) > 0)
  expect_true(sd(fit$test_predictions_sample_treat) > 0)
  
  # --- Check reproducibility ---
  fit2 <- CausalHorseForest(
    y = time,
    status = status,
    X_train_control = X_control,
    X_train_treat = X_treat,
    treatment_indicator_train = treatment,
    outcome_type = "right-censored",
    number_of_trees = 5,
    N_post = 10,
    N_burn = 5,
    store_posterior_sample = TRUE,
    verbose = FALSE,
    seed = 1
  )
  expect_equal(fit$train_predictions, fit2$train_predictions)
})


