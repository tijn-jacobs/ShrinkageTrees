test_that("calibrated horseshoe defaults are consistent across wrappers", {
  # The `k` defaults are literals in the wrapper signatures (so they render
  # readably in the usage sections), while the fallbacks in ShrinkageTrees()
  # and CausalShrinkageForest() use the constants. This test is what keeps the
  # two in sync.
  expect_equal(formals(HorseTrees)$k, ShrinkageTrees:::HORSESHOE_K_SINGLE)
  expect_equal(formals(CausalHorseForest)$k, ShrinkageTrees:::HORSESHOE_K_CAUSAL)
})

test_that("horseshoe scales default rather than erroring", {
  set.seed(1)
  n <- 40; p <- 5
  X <- matrix(runif(n * p), n, p)
  y <- as.numeric(2 * X[, 1] + rnorm(n, 0, 0.5))

  # Previously this combination raised
  #   "you must provide both local_hp and global_hp"
  expect_no_error(
    fit <- ShrinkageTrees(
      y = y, X_train = X, outcome_type = "continuous",
      number_of_trees = 5, prior_type = "horseshoe",
      N_post = 20, N_burn = 20, verbose = FALSE
    )
  )

  # The constructor nests the scales under $args, as param1 (local) and
  # param2 (global). See NewShrinkageTrees() in R/constructors.R.
  expected <- ShrinkageTrees:::HORSESHOE_K_SINGLE / sqrt(5)
  expect_equal(fit$args$param1, expected)
  expect_equal(fit$args$param2, expected)
})

test_that("supplying one scale leaves the other at its default", {
  set.seed(2)
  n <- 40; p <- 5
  X <- matrix(runif(n * p), n, p)
  y <- as.numeric(2 * X[, 1] + rnorm(n, 0, 0.5))

  fit <- ShrinkageTrees(
    y = y, X_train = X, outcome_type = "continuous",
    number_of_trees = 5, prior_type = "horseshoe",
    local_hp = 0.05,
    N_post = 20, N_burn = 20, verbose = FALSE
  )

  expect_equal(fit$args$param1, 0.05)
  expect_equal(fit$args$param2, ShrinkageTrees:::HORSESHOE_K_SINGLE / sqrt(5))
})

test_that("horseshoe_fw still requires both scales", {
  set.seed(3)
  n <- 40; p <- 5
  X <- matrix(runif(n * p), n, p)
  y <- as.numeric(2 * X[, 1] + rnorm(n, 0, 0.5))

  expect_error(
    ShrinkageTrees(
      y = y, X_train = X, outcome_type = "continuous",
      number_of_trees = 5, prior_type = "horseshoe_fw",
      N_post = 20, N_burn = 20, verbose = FALSE
    ),
    "horseshoe_fw"
  )
})

test_that("causal defaults use each forest's own tree count", {
  set.seed(4)
  n <- 60; p <- 5
  X <- matrix(runif(n * p), n, p)
  z <- rbinom(n, 1, 0.5)
  y <- as.numeric(X[, 1] + (z - 0.5) * 2 * X[, 2] + rnorm(n, 0, 0.5))

  expect_no_error(
    fit <- CausalShrinkageForest(
      y = y, X_train_control = X, X_train_treat = X,
      treatment_indicator_train = z,
      outcome_type = "continuous",
      number_of_trees_control = 5, number_of_trees_treat = 10,
      prior_type_control = "horseshoe", prior_type_treat = "horseshoe",
      N_post = 20, N_burn = 20, verbose = FALSE
    )
  )

  # The causal constructor nests per forest: $args$control and $args$treat.
  k <- ShrinkageTrees:::HORSESHOE_K_CAUSAL
  expect_equal(fit$args$control$param1, k / sqrt(5))
  expect_equal(fit$args$control$param2, k / sqrt(5))
  expect_equal(fit$args$treat$param1,   k / sqrt(10))
  expect_equal(fit$args$treat$param2,   k / sqrt(10))
})
