test_that("linear projection recovers the coefficients of exactly linear draws", {
  set.seed(1)
  n <- 40
  X <- cbind(a = runif(n), b = runif(n), c = runif(n))

  a0    <- c(1.0, 1.5, -0.5)
  betas <- rbind(c(2, -1, 0.5), c(1, 0, -2), c(-3, 2, 1))
  draws <- t(apply(cbind(a0, betas), 1, function(z) z[1] + X %*% z[-1]))

  pp <- posterior_projection(draws, X = X, family = "linear")

  expect_s3_class(pp, "PosteriorProjection")
  expect_equal(colnames(pp$coefficients), c("(Intercept)", "a", "b", "c"))
  expect_equal(unname(pp$coefficients[, 1]), a0, tolerance = 1e-4)
  expect_equal(unname(pp$coefficients[, -1]), betas, tolerance = 1e-4)

  # An exact linear summary leaves no residual, so the summary R^2 is 1.
  expect_equal(unname(pp$rho2), rep(1, 3), tolerance = 1e-6)
})

test_that("projection reports consistent dimensions", {
  set.seed(2)
  n <- 30
  S <- 8
  X <- cbind(a = rnorm(n), b = rnorm(n))
  draws <- matrix(rnorm(S * n), nrow = S)

  pp <- posterior_projection(draws, X = X, family = "linear")

  expect_equal(pp$n_obs, n)
  expect_equal(pp$n_draws, S)
  expect_equal(pp$p, 2)
  expect_equal(dim(pp$coefficients), c(S, 3))
  expect_length(pp$rho2, S)
  expect_length(pp$residuals, n)
  expect_equal(pp$covariates, c("a", "b"))
  expect_equal(nrow(pp$coef_summary), 3)
})

test_that("covariates selects a subset, by name or index", {
  set.seed(3)
  n <- 30
  X <- cbind(a = rnorm(n), b = rnorm(n), c = rnorm(n))
  draws <- matrix(rnorm(4 * n), nrow = 4)

  by_name <- posterior_projection(draws, X = X, covariates = c("a", "c"))
  by_idx  <- posterior_projection(draws, X = X, covariates = c(1, 3))

  expect_equal(by_name$covariates, c("a", "c"))
  expect_equal(by_name$p, 2)
  expect_equal(by_name$coefficients, by_idx$coefficients)
})

test_that("weights change the projection and are validated", {
  set.seed(4)
  n <- 30
  X <- cbind(a = rnorm(n))
  draws <- matrix(rnorm(3 * n), nrow = 3)

  flat <- posterior_projection(draws, X = X)
  wtd  <- posterior_projection(draws, X = X, weights = c(rep(3, 10), rep(1, 20)))

  expect_false(isTRUE(all.equal(flat$coefficients, wtd$coefficients)))
  expect_error(posterior_projection(draws, X = X, weights = rep(1, 5)),
               "weights must have length")
})

test_that("invalid input is rejected with an informative message", {
  set.seed(5)
  n <- 20
  X <- cbind(a = rnorm(n), b = rnorm(n))
  draws <- matrix(rnorm(3 * n), nrow = 3)

  expect_error(posterior_projection(draws, X = X, family = "additive",
                                    penalty = "lasso"),
               "applies only to family")

  bad <- draws
  bad[1, 1] <- NA
  expect_error(posterior_projection(bad, X = X), "non-finite")

  expect_error(posterior_projection(draws[, -1, drop = FALSE], X = X),
               "draws must be S x n")

  expect_error(posterior_projection(draws, X = X, covariates = "nope"),
               "not columns of the design")
})

test_that("additive family returns one effect curve per covariate", {
  set.seed(6)
  n <- 60
  X <- cbind(a = runif(n), b = runif(n))
  draws <- t(replicate(4, sin(3 * X[, 1]) + X[, 2] + rnorm(n, sd = 0.01)))

  pp <- posterior_projection(draws, X = X, family = "additive", df_spline = 3)

  expect_s3_class(pp, "PosteriorProjection")
  expect_equal(pp$family, "additive")
  expect_equal(pp$covariates, c("a", "b"))
  expect_true(all(is.finite(pp$rho2)))
  expect_gt(mean(pp$rho2), 0.5)
})

test_that("penalised projection summarises the posterior mean, without UQ", {
  skip_if_not_installed("glmnet")
  set.seed(7)
  n <- 50
  p <- 10
  X <- matrix(rnorm(n * p), n, p, dimnames = list(NULL, paste0("v", 1:p)))
  draws <- t(replicate(5, as.vector(2 * X[, 1] - X[, 2]) + rnorm(n, sd = 0.05)))

  pp <- posterior_projection(draws, X = X, family = "linear",
                             penalty = "lasso", lambda = "min")

  expect_equal(pp$penalty, "lasso")
  expect_true(!is.null(pp$structure$cv))

  # One fit on the Bayes estimate: one row of coefficients, one rho^2, and no
  # interval columns anywhere.
  expect_false(pp$uq)
  expect_equal(nrow(pp$coefficients), 1L)
  expect_length(pp$rho2, 1L)
  expect_true(all(is.na(pp$coef_summary$lower)))
  expect_true(all(is.na(pp$coef_summary$upper)))
  expect_true(all(is.na(pp$coef_summary$sd)))
  expect_null(pp$structure$nonzero)

  # Only selected covariates are reported, and the signal is among them.
  expect_true(all(pp$coefficients[1, -1] != 0))
  expect_true(all(c("v1", "v2") %in% pp$covariates))

  # The projection is of the posterior mean, so it does not depend on the
  # order of the draws or on anything else about their spread.
  pp2 <- posterior_projection(draws[c(3, 1, 5, 2, 4), ], X = X,
                              family = "linear", penalty = "lasso",
                              lambda = pp$structure$lambda)
  pp1 <- posterior_projection(draws, X = X, family = "linear",
                              penalty = "lasso",
                              lambda = pp$structure$lambda)
  expect_equal(pp2$coefficients, pp1$coefficients)
})

test_that("unpenalised projection keeps its uncertainty", {
  set.seed(17)
  n <- 40
  X <- cbind(a = rnorm(n), b = rnorm(n))
  draws <- t(replicate(20, as.vector(X %*% c(1, -2)) + rnorm(n, sd = 0.2)))

  pp <- posterior_projection(draws, X = X, family = "linear")

  expect_true(pp$uq)
  expect_equal(nrow(pp$coefficients), 20L)
  expect_true(all(is.finite(pp$coef_summary$lower)))
  expect_true(all(pp$coef_summary$lower <= pp$coef_summary$upper))
})

test_that("tree family partitions into at most max_leaves groups", {
  skip_if_not_installed("rpart")
  set.seed(8)
  n <- 80
  X <- cbind(a = runif(n), b = runif(n))
  draws <- t(replicate(4, ifelse(X[, 1] > 0.5, 2, -2) + rnorm(n, sd = 0.05)))

  pp <- posterior_projection(draws, X = X, family = "tree", max_leaves = 4)

  expect_equal(pp$family, "tree")
  expect_lte(ncol(pp$coefficients), 4)
  expect_gt(mean(pp$rho2), 0.5)
})

test_that("projection dispatches on a fitted ShrinkageTrees object", {
  set.seed(9)
  n <- 40
  p <- 3
  X <- matrix(runif(n * p), n, p, dimnames = list(NULL, c("a", "b", "c")))
  y <- X[, 1] + rnorm(n, sd = 0.1)

  fit <- ShrinkageTrees(y = y, X_train = X, outcome_type = "continuous",
                        number_of_trees = 5, prior_type = "horseshoe",
                        N_post = 10, N_burn = 5,
                        store_posterior_sample = TRUE, verbose = FALSE)

  pp <- posterior_projection(fit, family = "linear")

  expect_s3_class(pp, "PosteriorProjection")
  expect_equal(pp$n_obs, n)
  expect_equal(pp$n_draws, nrow(fit$train_predictions_sample))
  expect_equal(pp$covariates, c("a", "b", "c"))
  expect_equal(pp$scale, "response")
  expect_equal(pp$target, "f")
})

test_that("projection needs a stored posterior sample", {
  set.seed(10)
  n <- 30
  X <- matrix(runif(n * 2), n, 2)
  y <- X[, 1] + rnorm(n, sd = 0.1)

  fit <- ShrinkageTrees(y = y, X_train = X, outcome_type = "continuous",
                        number_of_trees = 5, prior_type = "horseshoe",
                        N_post = 10, N_burn = 5,
                        store_posterior_sample = FALSE, verbose = FALSE)

  expect_error(posterior_projection(fit), "store_posterior_sample")
})

test_that("print returns its argument invisibly", {
  set.seed(11)
  n <- 30
  X <- cbind(a = rnorm(n), b = rnorm(n))
  draws <- matrix(rnorm(4 * n), nrow = 4)
  pp <- posterior_projection(draws, X = X)

  out <- withVisible(print(pp))
  expect_false(out$visible)
  expect_identical(out$value, pp)
})

test_that("plot rejects types the projection cannot supply", {
  skip_if_not_installed("ggplot2")
  set.seed(12)
  n <- 30
  X <- cbind(a = rnorm(n), b = rnorm(n))
  draws <- matrix(rnorm(4 * n), nrow = 4)
  pp <- posterior_projection(draws, X = X, family = "linear")

  expect_error(plot(pp, type = "path"), "cross-validated penalty")
  expect_error(plot(pp, type = "tree"), "family = 'tree'")
  expect_error(plot(pp, type = "effects"), "family = 'additive'")

  expect_s3_class(plot(pp, type = "coefficients"), "ggplot")
  expect_s3_class(plot(pp, type = "rho2"), "ggplot")
})

test_that("plot refuses an rho2 density when there is no posterior", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("glmnet")
  set.seed(13)
  n <- 40
  p <- 8
  X <- matrix(rnorm(n * p), n, p, dimnames = list(NULL, paste0("v", 1:p)))
  draws <- t(replicate(5, as.vector(X[, 1]) + rnorm(n, sd = 0.05)))

  pp <- posterior_projection(draws, X = X, family = "linear",
                             penalty = "lasso")

  expect_error(plot(pp, type = "rho2"), "projection per draw")
  expect_s3_class(plot(pp, type = "coefficients"), "ggplot")
  expect_s3_class(plot(pp, type = "path"), "ggplot")
})
