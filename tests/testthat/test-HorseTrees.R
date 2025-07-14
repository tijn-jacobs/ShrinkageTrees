test_that("HorseTrees for continuous outcome", {
  
  # Generate data
  X <- matrix(runif(50 * 3), ncol = 3)
  X_test <- matrix(runif(50 * 3), ncol = 3)
  y <- X[, 1] + rnorm(50)
  
  # Fit the model
  fit <- HorseTrees(y = y, 
                    X_train = X,
                    X_test = X_test,
                    outcome_type = "continuous",
                    number_of_trees = 5,
                    N_post = 10,
                    N_burn = 5,
                    store_posterior_sample = TRUE,
                    verbose = FALSE,
                    seed = 1)
  
  # Basic checks
  expect_type(fit, "list")
  expect_true("train_predictions" %in% names(fit))
  expect_true("test_predictions" %in% names(fit))
  expect_length(fit$train_predictions, 50)
  expect_length(fit$sigma, 10 + 5)
  
  # Numerical checks
  expect_false(any(is.na(fit$train_predictions)))
  expect_false(any(is.nan(fit$train_predictions)))
  expect_true(all(is.finite(fit$train_predictions)))
  
  # Check posterior samples
  expect_true("train_predictions_sample" %in% names(fit))
  expect_true("test_predictions_sample" %in% names(fit))
  expect_equal(dim(fit$train_predictions_sample), c(10, 50))
  expect_equal(dim(fit$test_predictions_sample), c(10, 50))
  
  # Check sigma positivity
  expect_true(all(fit$sigma > 0))
  
  # Check predictions are not all zero
  expect_true(sd(fit$train_predictions) > 0)
  expect_true(sd(fit$train_predictions_sample) > 0)
  expect_true(sd(fit$test_predictions) > 0)
  expect_true(sd(fit$test_predictions_sample) > 0)
  
  # Check reproducibility
  fit2 <- HorseTrees(y = y, X_train = X, outcome_type = "continuous", number_of_trees = 5, 
                     N_post = 10, N_burn = 5, store_posterior_sample = TRUE, 
                     verbose = FALSE, seed = 1)
  expect_equal(fit$train_predictions, fit2$train_predictions)
})


test_that("HorseTrees for binary outcome", {
  
  # Generate data
  X <- matrix(rnorm(50 * 3), ncol = 3)
  X_test <- matrix(runif(50 * 3), ncol = 3)
  y <- ifelse(X[, 1] + rnorm(50) > 0, 1, 0)
  
  # Fit the model
  fit <- HorseTrees(y = y, 
                    X_train = X,
                    X_test = X_test,
                    outcome_type = "binary",
                    number_of_trees = 5,
                    N_post = 10,
                    N_burn = 5,
                    store_posterior_sample = TRUE,
                    verbose = FALSE,
                    seed = 1)
  
  # Basic checks
  expect_type(fit, "list")
  expect_true("train_predictions" %in% names(fit))
  expect_length(fit$train_predictions, 50)

  # Numerical checks
  expect_false(any(is.na(fit$train_predictions)))
  expect_false(any(is.nan(fit$train_predictions)))
  expect_true(all(is.finite(fit$train_predictions)))
  
  # Check posterior samples
  expect_true("train_predictions_sample" %in% names(fit))
  expect_true("test_predictions_sample" %in% names(fit))
  expect_true("train_probabilities_sample" %in% names(fit))
  expect_true("test_probabilities_sample" %in% names(fit))
  expect_equal(dim(fit$train_predictions_sample), c(10, 50))
  expect_equal(dim(fit$test_predictions_sample), c(10, 50))
  expect_equal(dim(fit$train_probabilities_sample), c(10, 50))
  expect_equal(dim(fit$test_probabilities_sample), c(10, 50))
  
  # Check estimated probabilitiesHmm ik are within [0,1]
  expect_true(all(fit$train_probabilities >= 0 & fit$train_probabilities <= 1))
  expect_true(all(fit$test_probabilities >= 0 & fit$test_probabilities <= 1))
  
  # Check predictions are not all constant
  expect_true(sd(fit$train_predictions) > 0)
  expect_true(sd(fit$train_predictions_sample) > 0)
  expect_true(sd(fit$test_predictions) > 0)
  expect_true(sd(fit$test_predictions_sample) > 0)
  
  # Check reproducibility
  fit2 <- HorseTrees(y = y, X_train = X, outcome_type = "binary", number_of_trees = 5, 
                     N_post = 10, N_burn = 5, store_posterior_sample = TRUE, 
                     verbose = FALSE, seed = 1)
  expect_equal(fit$train_predictions, fit2$train_predictions)
})

test_that("HorseTrees for right-censored survival outcome", {
  
  # Generate covariates
  X <- matrix(rnorm(50 * 3), ncol = 3)
  X_test <- matrix(runif(50 * 3), ncol = 3)
  
  # Generate survival times (exponential baseline hazard depending on X[,1])
  linpred <- X[, 1]
  true_time <- rexp(50, rate = exp(linpred))
  
  # Generate censoring times
  censor_time <- rexp(50, rate = 0.5)
  
  # Observed time and status indicator
  time_obs <- pmin(true_time, censor_time)
  status <- as.numeric(true_time <= censor_time)
  
  # Fit the model
  fit <- HorseTrees(y = time_obs, 
                    X_train = X,
                    X_test = X_test,
                    outcome_type = "right-censored",
                    status = status,
                    number_of_trees = 5,
                    N_post = 10,
                    N_burn = 5,
                    store_posterior_sample = TRUE,
                    verbose = FALSE,
                    seed = 1)
  
  # Basic checks
  expect_type(fit, "list")
  expect_true("train_predictions" %in% names(fit))
  expect_true("test_predictions" %in% names(fit))
  expect_true("test_predictions" %in% names(fit))
  expect_length(fit$train_predictions, 50)
  expect_length(fit$sigma, 10 + 5)
  
  # Numerical checks
  expect_false(any(is.na(fit$train_predictions)))
  expect_false(any(is.nan(fit$train_predictions)))
  expect_true(all(is.finite(fit$train_predictions)))
  
  # Check posterior samples
  expect_true("train_predictions_sample" %in% names(fit))
  expect_true("test_predictions_sample" %in% names(fit))
  expect_equal(dim(fit$train_predictions_sample), c(10, 50))
  expect_equal(dim(fit$test_predictions_sample), c(10, 50))
  
  # Check sigma positivity
  expect_true(all(fit$sigma > 0))
  
  # Check predictions are not all constant
  expect_true(sd(fit$train_predictions) > 0)
  expect_true(sd(fit$train_predictions_sample) > 0)
  expect_true(sd(fit$test_predictions) > 0)
  expect_true(sd(fit$test_predictions_sample) > 0)
  
  # Check reproducibility
  fit2 <- HorseTrees(y = time_obs, X_train = X, outcome_type = "right-censored", 
                     status = status, number_of_trees = 5, N_post = 10, 
                     N_burn = 5, store_posterior_sample = TRUE, verbose = FALSE, 
                     seed = 1)
  expect_equal(fit$train_predictions, fit2$train_predictions)
})
