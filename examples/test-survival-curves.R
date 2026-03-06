# Quick test: survival curve monotonicity
library(ShrinkageTrees)

set.seed(1)
n <- 80
X <- matrix(rnorm(n * 5), ncol = 5)
time <- exp(1 + 0.5 * X[, 1] + rnorm(n, 0, 0.5))
status <- rbinom(n, 1, 0.7)

# --- Training curves ---
fit <- SurvivalBART(
  time = time, status = status, X_train = X,
  number_of_trees = 200, N_post = 5000, N_burn = 2500,
  verbose = FALSE
)

p1 <- plot(fit, type = "survival", xlim=c(0,15))
print(p1)

p2 <- plot(fit, type = "survival", obs = c(1, 10, 20), xlim=c(0,15))
print(p2)

# --- Prediction curves ---
X_new <- matrix(rnorm(20 * 5), ncol = 5)
pred <- predict(fit, newdata = X_new)

p3 <- plot(pred, type = "survival", xlim=c(0,15))
print(p3)

p4 <- plot(pred, type = "survival", obs = c(1, 10, 20), xlim=c(0,15))
print(p4)
