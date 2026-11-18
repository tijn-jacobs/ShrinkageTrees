# Simulated data for tests/testthat/test-refactor-equivalence.R.
#
# One data set per seed: a design matrix, a test matrix, a treatment
# indicator, and continuous, binary, causal, right-censored and
# interval-censored outcomes.

refactor_data <- function(seed, n = 40, p = 4) {
  set.seed(seed)
  X      <- matrix(runif(n * p), n, p)
  X_test <- matrix(runif(10 * p), 10, p)
  trt    <- rbinom(n, 1, 0.5)
  mu     <- 2 * X[, 1] - X[, 2]

  y_cont <- as.numeric(mu + rnorm(n, 0, 0.5))
  y_bin  <- as.numeric(mu + rnorm(n) > 0.5)

  # Right-censored survival times
  t_event <- exp(0.5 * mu + trt * 0.4 + rnorm(n, 0, 0.4))
  t_cens  <- exp(rnorm(n, 0.8, 0.5))
  time    <- pmin(t_event, t_cens)
  status  <- as.numeric(t_event <= t_cens)

  # Interval-censored version: a third exact, a third right-censored,
  # a third observed in a finite interval
  kind  <- rep(1:3, length.out = n)
  left  <- ifelse(kind == 3, t_event * 0.8, t_event)
  right <- ifelse(kind == 1, t_event, ifelse(kind == 2, Inf, t_event * 1.3))

  list(X = X, X_test = X_test, trt = trt, y_cont = y_cont, y_bin = y_bin,
       y_causal = as.numeric(y_cont + trt * (1 + X[, 3])),
       time = time, status = status, left = left, right = right,
       propensity = rep(0.5, n), n = n, p = p)
}
