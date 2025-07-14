test_that("HorseTrees works for continuous outcome", {
  set.seed(1)
  X <- matrix(rnorm(50 * 3), ncol = 3)
  y <- X[, 1] + rnorm(50)
  
  fit <- HorseTrees(y, X,
                    outcome_type = "continuous",
                    number_of_trees = 5,
                    N_post = 10,
                    N_burn = 5,
                    store_posterior_sample = TRUE,
                    verbose = FALSE)
  
  expect_type(fit, "list")
  expect_true("train_predictions" %in% names(fit))
  expect_length(fit$train_predictions, 50)
  expect_length(fit$sigma, 10 + 5)
})

test_that("HorseTrees works for binary outcome", {
  set.seed(2)
  X <- matrix(rnorm(50 * 3), ncol = 3)
  y <- ifelse(X[, 1] + rnorm(50) > 0, 1, 0)
  
  fit <- HorseTrees(y, X,
                    outcome_type = "binary",
                    number_of_trees = 5,
                    N_post = 10,
                    N_burn = 5,
                    store_posterior_sample = FALSE,
                    verbose = FALSE)
  
  expect_type(fit, "list")
  expect_true("train_predictions" %in% names(fit))
  expect_length(fit$train_predictions, 50)
})

test_that("HorseTrees works for right-censored outcome", {
  set.seed(3)
  X <- matrix(rnorm(50 * 3), ncol = 3)
  time <- rexp(50, rate = 0.1)
  status <- rbinom(50, 1, 0.7)
  
  fit <- HorseTrees(time, X,
                    status = status,
                    outcome_type = "right-censored",
                    timescale = "time",
                    number_of_trees = 5,
                    N_post = 10,
                    N_burn = 5,
                    store_posterior_sample = FALSE,
                    verbose = FALSE)
  
  expect_type(fit, "list")
  expect_true("train_predictions" %in% names(fit))
  expect_length(fit$train_predictions, 50)
})
