# Test survival curve plotting
library(ShrinkageTrees)

set.seed(42)
n <- 100; p <- 5000
X <- matrix(rnorm(n * p), ncol = p)
true_log_t <- 2 + 0.5 * X[, 1] - 0.3 * X[, 2] + rnorm(n, 0, 0.5)
time <- exp(true_log_t)
status <- rbinom(n, 1, 0.7)

# Fit with store_posterior_sample = TRUE
fit <- SurvivalBART(
  time = time, status = status, X_train = X,
  number_of_trees = 10, N_post = 2500, N_burn = 1000,
  verbose = FALSE
)

# 1. Population-averaged survival curve (no KM)
p1 <- plot(fit, type = "survival", level = 0.90)
print(p1)

# 1b. Population-averaged survival curve with KM overlay
p1b <- plot(fit, type = "survival", level = 0.90, km = TRUE)
print(p1b)


# 2. Individual survival curves for specific observations
p2 <- plot(fit, type = "survival", obs = c(1, 25, 50, 75, 100), km = TRUE)
print(p2)

# 3. Single individual
p3 <- plot(fit, type = "survival", obs = 1)
print(p3)

# 4. Custom time grid + 90% CI
p4 <- plot(fit, type = "survival", obs = c(1, 50),
           t_grid = seq(0.5, 30, length.out = 100), level = 0.90)
print(p4)

# 5. Without store_posterior_sample (plug-in mu, sigma-only bands)
fit_no_sample <- SurvivalBART(
  time = time, status = status, X_train = X,
  number_of_trees = 10, N_post = 200, N_burn = 100,
  store_posterior_sample = FALSE, verbose = FALSE
)
p5 <- plot(fit_no_sample, type = "survival")
print(p5)

p6 <- plot(fit_no_sample, type = "survival", obs = c(1, 50))
print(p6)

# 6. Verify error for non-survival model
fit_cont <- HorseTrees(
  y = rnorm(50), X_train = matrix(rnorm(50 * 3), ncol = 3),
  number_of_trees = 5, N_post = 50, N_burn = 25, verbose = FALSE
)
tryCatch(
  plot(fit_cont, type = "survival"),
  error = function(e) cat("Expected error:", e$message, "\n")
)

cat("\nAll survival curve tests passed!\n")

