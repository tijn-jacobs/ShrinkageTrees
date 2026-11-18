# Tests for the 2.1.1 clean-up.
#
# They pin the relations the clean-up relies on: HorseTrees() and
# CausalHorseForest() are special cases of the general fitters, and the
# shared helpers behave as the inlined code did. Each comparison runs both
# sides under the same seed in the same session, so it holds for any random
# number generator.

test_that("HorseTrees equals ShrinkageTrees with the horseshoe prior", {
  d  <- refactor_data(201)
  k  <- 0.8
  m  <- 5
  hp <- k / sqrt(m)

  for (type in c("continuous", "right-censored")) {
    y      <- if (type == "continuous") d$y_cont else d$time
    status <- if (type == "continuous") NULL else d$status

    set.seed(1)
    fit_h <- HorseTrees(
      y = y, status = status, X_train = d$X, X_test = d$X_test,
      outcome_type = type, number_of_trees = m, k = k,
      N_post = 20, N_burn = 10, verbose = FALSE
    )
    set.seed(1)
    fit_s <- ShrinkageTrees(
      y = y, status = status, X_train = d$X, X_test = d$X_test,
      outcome_type = type, number_of_trees = m,
      prior_type = "horseshoe", local_hp = hp, global_hp = hp,
      N_post = 20, N_burn = 10, verbose = FALSE
    )

    expect_identical(fit_h$train_predictions, fit_s$train_predictions)
    expect_identical(fit_h$test_predictions,  fit_s$test_predictions)
    expect_identical(fit_h$sigma,             fit_s$sigma)
    expect_identical(fit_h$train_predictions_sample,
                     fit_s$train_predictions_sample)
    expect_identical(fit_h$args, fit_s$args)
  }

  # The wrapper reports its own call, not the inner one
  expect_identical(as.character(fit_h$call[[1]]), "HorseTrees")
  expect_s3_class(fit_h, "ShrinkageTrees")
})

test_that("CausalHorseForest equals CausalShrinkageForest with horseshoe priors", {
  d  <- refactor_data(202)
  k  <- 1.2
  m  <- 5
  hp <- k / sqrt(m)

  set.seed(1)
  fit_h <- CausalHorseForest(
    y = d$y_causal, X_train_control = d$X, X_train_treat = d$X,
    treatment_indicator_train = d$trt, outcome_type = "continuous",
    number_of_trees = m, k = k, power = 2.5, base = 0.9,
    N_post = 20, N_burn = 10, store_posterior_sample = TRUE, verbose = FALSE
  )
  set.seed(1)
  fit_s <- CausalShrinkageForest(
    y = d$y_causal, X_train_control = d$X, X_train_treat = d$X,
    treatment_indicator_train = d$trt, outcome_type = "continuous",
    number_of_trees_control = m, number_of_trees_treat = m,
    prior_type_control = "horseshoe", prior_type_treat = "horseshoe",
    local_hp_control = hp, local_hp_treat = hp,
    global_hp_control = hp, global_hp_treat = hp,
    power_control = 2.5, power_treat = 2.5,
    base_control = 0.9, base_treat = 0.9,
    # CausalHorseForest() defaults; CausalShrinkageForest() defaults to 0.5
    p_grow = 0.4, p_prune = 0.4,
    N_post = 20, N_burn = 10, store_posterior_sample = TRUE, verbose = FALSE
  )

  expect_identical(fit_h$train_predictions,       fit_s$train_predictions)
  expect_identical(fit_h$train_predictions_treat, fit_s$train_predictions_treat)
  expect_identical(fit_h$sigma,                   fit_s$sigma)
  expect_identical(fit_h$train_predictions_sample_treat,
                   fit_s$train_predictions_sample_treat)
  expect_identical(fit_h$args, fit_s$args)

  expect_identical(as.character(fit_h$call[[1]]), "CausalHorseForest")
  expect_s3_class(fit_h, "CausalShrinkageForest")
})

test_that("multi-chain fits report the call of the function the user called", {
  d <- refactor_data(203)
  fit <- HorseTrees(
    y = d$y_cont, X_train = d$X, outcome_type = "continuous",
    number_of_trees = 5, N_post = 10, N_burn = 5, n_chains = 2,
    verbose = FALSE
  )
  expect_identical(as.character(fit$call[[1]]), "HorseTrees")
  expect_equal(fit$mcmc$n_chains, 2)
  expect_length(fit$sigma, 20)
  expect_equal(dim(fit$train_predictions_sample), c(20, d$n))
})

test_that(".interval_censoring encodes exact, right- and interval-censored rows", {
  left  <- c(1, 2, 3)
  right <- c(1, Inf, 5)
  ic <- .interval_censoring(left, right)

  expect_identical(ic$status,       c(1L, 0L, 0L))
  expect_identical(ic$ic_indicator, c(0L, 0L, 1L))
  expect_equal(ic$y, c(1, 2, 4))
})

test_that("censored_info returns only the centring mean and scale", {
  skip_if_not_installed("survival")
  d <- refactor_data(204)
  info <- censored_info(log(d$time), d$status)

  expect_named(info, c("mu", "sd"))
  ref <- survival::survreg(survival::Surv(log(d$time), d$status) ~ 1,
                           dist = "gaussian")
  expect_equal(info$mu, unname(ref$coefficients[1]))
  expect_equal(info$sd, ref$scale)
})

test_that("continuous causal fits no longer carry empty sample elements", {
  d <- refactor_data(205)
  fit <- CausalShrinkageForest(
    y = d$y_causal, X_train_control = d$X, X_train_treat = d$X,
    treatment_indicator_train = d$trt, outcome_type = "continuous",
    number_of_trees_control = 5, number_of_trees_treat = 5,
    N_post = 10, N_burn = 5, store_posterior_sample = TRUE, verbose = FALSE
  )
  expect_false("train_predictions_sample" %in% names(fit))
  expect_false("test_predictions_sample"  %in% names(fit))
  expect_equal(dim(fit$train_predictions_sample_control), c(10, d$n))
})
