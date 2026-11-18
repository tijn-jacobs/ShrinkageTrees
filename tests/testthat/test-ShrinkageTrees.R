test_that("ShrinkageTrees (horseshoe) works for continuous outcome", {
  
  # Generate data
  X <- matrix(runif(50 * 3), ncol = 3)
  X_test <- matrix(runif(50 * 3), ncol = 3)
  y <- X[, 1] + rnorm(50)
  
  # Fit the model
  set.seed(1)
  fit <- ShrinkageTrees(y = y, 
                        X_train = X,
                        X_test = X_test,
                        outcome_type = "continuous",
                        number_of_trees = 50,
                        prior_type = "horseshoe",
                        local_hp = 0.1/sqrt(50),
                        global_hp = 0.1/sqrt(50),
                        N_post = 10,
                        N_burn = 5,
                        store_posterior_sample = TRUE,
                        verbose = FALSE)
  
  # Basic checks
  expect_type(fit, "list")
  expect_true("train_predictions" %in% names(fit))
  expect_true("test_predictions" %in% names(fit))
  expect_length(fit$train_predictions, 50)
  expect_length(fit$sigma, 10)
  
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
  set.seed(1)
  fit2 <- ShrinkageTrees(y = y, 
                         X_train = X,
                         X_test = X_test,
                         outcome_type = "continuous",
                         number_of_trees = 50,
                         prior_type = "horseshoe",
                         local_hp = 0.1/sqrt(50),
                         global_hp = 0.1/sqrt(50),
                         N_post = 10,
                         N_burn = 5,
                         store_posterior_sample = TRUE,
                         verbose = FALSE)
  expect_equal(fit$train_predictions, fit2$train_predictions)
})


test_that("ShrinkageTrees (horseshoe) works for binary outcome", {
  
  # Generate data
  X <- matrix(rnorm(50 * 3), ncol = 3)
  X_test <- matrix(runif(50 * 3), ncol = 3)
  y <- ifelse(X[, 1] + rnorm(50) > 0, 1, 0)
  
  # Fit the model
  set.seed(1)
  fit <- ShrinkageTrees(y = y, 
                        X_train = X,
                        X_test = X_test,
                        outcome_type = "binary",
                        number_of_trees = 50,
                        prior_type = "horseshoe",
                        local_hp = 0.1/sqrt(50),
                        global_hp = 0.1/sqrt(50),
                        N_post = 10,
                        N_burn = 5,
                        store_posterior_sample = TRUE,
                        verbose = FALSE)
    
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
  
  # Check estimated probabilities are within [0,1]
  expect_true(all(fit$train_probabilities >= 0 & fit$train_probabilities <= 1))
  expect_true(all(fit$test_probabilities >= 0 & fit$test_probabilities <= 1))
  
  # Check predictions are not all constant
  expect_true(sd(fit$train_predictions) > 0)
  expect_true(sd(fit$train_predictions_sample) > 0)
  expect_true(sd(fit$test_predictions) > 0)
  expect_true(sd(fit$test_predictions_sample) > 0)
  
  # Check reproducibility
  set.seed(1)
  fit2 <- ShrinkageTrees(y = y, 
                         X_train = X,
                         X_test = X_test,
                         outcome_type = "binary",
                         number_of_trees = 50,
                         prior_type = "horseshoe",
                         local_hp = 0.1/sqrt(50),
                         global_hp = 0.1/sqrt(50),
                         N_post = 10,
                         N_burn = 5,
                         store_posterior_sample = TRUE,
                         verbose = FALSE)
  expect_equal(fit$train_predictions, fit2$train_predictions)
})

test_that("ShrinkageTrees (horseshoe) works for right-censored survival outcome", {
  
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
  set.seed(1)
  fit <- ShrinkageTrees(y = time_obs, 
                        X_train = X,
                        X_test = X_test,
                        outcome_type = "right-censored",
                        status = status,
                        number_of_trees = 50,
                        prior_type = "horseshoe",
                        local_hp = 0.1/sqrt(50),
                        global_hp = 0.1/sqrt(50),
                        N_post = 10,
                        N_burn = 5,
                        store_posterior_sample = TRUE,
                        verbose = FALSE)
  
  # Basic checks
  expect_type(fit, "list")
  expect_true("train_predictions" %in% names(fit))
  expect_true("test_predictions" %in% names(fit))
  expect_true("test_predictions" %in% names(fit))
  expect_length(fit$train_predictions, 50)
  expect_length(fit$sigma, 10)
  
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
  set.seed(1)
  fit2 <- ShrinkageTrees(y = time_obs, 
                         X_train = X,
                         X_test = X_test,
                         outcome_type = "right-censored",
                         status = status,
                         number_of_trees = 50,
                         prior_type = "horseshoe",
                         local_hp = 0.1/sqrt(50),
                         global_hp = 0.1/sqrt(50),
                         N_post = 10,
                         N_burn = 5,
                         store_posterior_sample = TRUE,
                         verbose = FALSE)
  expect_equal(fit$train_predictions, fit2$train_predictions)
})


test_that("ShrinkageTrees (half-cauchy) works for continuous outcome", {
  
  # Generate data
  X <- matrix(runif(50 * 3), ncol = 3)
  X_test <- matrix(runif(50 * 3), ncol = 3)
  y <- X[, 1] + rnorm(50)
  
  # Fit the model
  set.seed(1)
  fit <- ShrinkageTrees(y = y, 
                        X_train = X,
                        X_test = X_test,
                        outcome_type = "continuous",
                        number_of_trees = 50,
                        prior_type = "half-cauchy",
                        local_hp = 0.1/sqrt(50),
                        N_post = 10,
                        N_burn = 5,
                        store_posterior_sample = TRUE,
                        verbose = FALSE)
  
  # Basic checks
  expect_type(fit, "list")
  expect_true("train_predictions" %in% names(fit))
  expect_true("test_predictions" %in% names(fit))
  expect_length(fit$train_predictions, 50)
  expect_length(fit$sigma, 10)
  
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
  set.seed(1)
  fit2 <- ShrinkageTrees(y = y, 
                         X_train = X,
                         X_test = X_test,
                         outcome_type = "continuous",
                         number_of_trees = 50,
                         prior_type = "half-cauchy",
                         local_hp = 0.1/sqrt(50),
                         N_post = 10,
                         N_burn = 5,
                         store_posterior_sample = TRUE,
                         verbose = FALSE)
  expect_equal(fit$train_predictions, fit2$train_predictions)
})

test_that("ShrinkageTrees errors on invalid outcome_type", {
  X <- matrix(rnorm(20 * 2), ncol = 2)
  y <- rnorm(20)
  
  expect_error(
    ShrinkageTrees(
      y = y,
      X_train = X,
      outcome_type = "invalid"
    )
  )
})

test_that("X_test defaults correctly", {
  X <- matrix(rnorm(40), ncol = 2)
  y <- rnorm(20)

  fit <- ShrinkageTrees(
    y = y,
    X_train = X,
    outcome_type = "continuous",
    number_of_trees = 5,
    prior_type = "standard",
    local_hp = 0.1,
    N_post = 5,
    N_burn = 2,
    verbose = FALSE
  )

  expect_length(fit$test_predictions, 1)
})


# ── S3 class, print, summary, predict ────────────────────────────────────────

test_that("ShrinkageTrees S3 methods work", {
  X      <- matrix(runif(50 * 3), ncol = 3)
  X_test <- matrix(runif(10 * 3), ncol = 3)
  y      <- X[, 1] + rnorm(50)

  set.seed(1)
  fit <- ShrinkageTrees(
    y = y, X_train = X,
    outcome_type = "continuous",
    number_of_trees = 5,
    prior_type = "horseshoe",
    local_hp = 0.1 / sqrt(5), global_hp = 0.1 / sqrt(5),
    N_post = 10, N_burn = 5,
    store_posterior_sample = TRUE,
    verbose = FALSE
  )

  expect_s3_class(fit, "ShrinkageTrees")
  expect_no_error(capture.output(print(fit)))
  expect_no_error(smry <- summary(fit))
  expect_type(smry, "list")
  expect_true(!is.null(smry$sigma))

  # predict on new data
  expect_no_error(pred <- predict(fit, newdata = X_test))
  expect_s3_class(pred, "ShrinkageTreesPrediction")
  expect_length(pred$mean, nrow(X_test))
  expect_length(pred$lower, nrow(X_test))
  expect_length(pred$upper, nrow(X_test))
  expect_true(all(pred$lower <= pred$mean))
  expect_true(all(pred$mean  <= pred$upper))

  expect_no_error(capture.output(print(pred)))
  expect_no_error(capture.output(summary(pred)))
})


# ── Multi-chain (n_chains) ────────────────────────────────────────────────────

test_that("ShrinkageTrees n_chains > 1 pools chains correctly", {
  X <- matrix(runif(50 * 3), ncol = 3)
  y <- X[, 1] + rnorm(50)

  set.seed(1)
  fit <- ShrinkageTrees(
    y = y, X_train = X,
    outcome_type = "continuous",
    number_of_trees = 5,
    prior_type = "horseshoe",
    local_hp = 0.1 / sqrt(5), global_hp = 0.1 / sqrt(5),
    N_post = 10, N_burn = 5,
    n_chains = 2,
    verbose = FALSE
  )

  expect_s3_class(fit, "ShrinkageTrees")
  expect_length(fit$sigma, 20)
  expect_equal(fit$mcmc$n_chains, 2)
  expect_length(fit$chains$acceptance_ratios, 2)
})


# ── summary() reports test predictions only when X_test was given ─────────────

test_that("summary reports test predictions only when X_test was given", {
  set.seed(1)
  n <- 30; p <- 3
  X      <- matrix(runif(n * p), ncol = p)
  X_test <- matrix(runif(5 * p), ncol = p)
  y      <- X[, 1] + rnorm(n)

  fit_no_test <- HorseTrees(y = y, X_train = X, number_of_trees = 5,
                            N_post = 10, N_burn = 5, verbose = FALSE)
  expect_false(fit_no_test$data_info$test_provided)
  expect_null(summary(fit_no_test)$predictions$test)
  expect_false(any(grepl("Test:", capture.output(print(summary(fit_no_test))))))

  fit_test <- ShrinkageTrees(y = y, X_train = X, X_test = X_test,
                             prior_type = "horseshoe", number_of_trees = 5,
                             N_post = 10, N_burn = 5, verbose = FALSE)
  expect_true(fit_test$data_info$test_provided)
  expect_false(is.null(summary(fit_test)$predictions$test))
  expect_true(any(grepl("Test:", capture.output(print(summary(fit_test))))))

  # Multi-chain fits keep the flag of the per-chain fits
  fit_chains <- HorseTrees(y = y, X_train = X, number_of_trees = 5,
                           N_post = 10, N_burn = 5, n_chains = 2,
                           verbose = FALSE)
  expect_false(fit_chains$data_info$test_provided)
})
