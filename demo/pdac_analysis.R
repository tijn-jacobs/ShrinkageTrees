## Demo: Analysis of the TCGA PAAD dataset (pdac)
## ShrinkageTrees package

# Load package
library(ShrinkageTrees)

# Load the data
data("pdac")

# Retrieve the data
time <- pdac$time
status <- pdac$status
treatment <- pdac$treatment
covariates <- pdac[, !(colnames(pdac) %in% c("time", "status", "treatment"))]

# Estimate the propensity scores
propensity_fit <- ShrinkageTrees::HorseTrees(
  y = treatment,
  X_train = covariates,
  outcome_type = "binary",
  k = 0.1,
  N_post = 5000,
  N_burn = 5000
)

# Retrieve estimated propensity scores
propensity <- pnorm(propensity_fit$train_predictions)

# Define p0 and p1
p0 <- propensity[treatment == 0]
p1 <- propensity[treatment == 1]

# Define colors
col_treated <- rgb(0, 0.5, 0, 0.5)  # green
col_control <- rgb(1, 0.5, 0, 0.5)  # orange

# Plot histograms
hist(p0,
     breaks = 10,
     xlim = c(0.33, 0.63),
     col = col_control,
     xlab = "Propensity score",
     ylab = "Frequency",
     main = "Propensity score overlap",
)
hist(p1,
     breaks = 10,
     xlim = range(propensity),
     col = col_treated,
     add = TRUE
)
legend("topright",
       legend = c("Control", "Treated"),
       fill = c(col_control, col_treated)
)

# Placeholder for your analysis
cat("\nAnalysis section TBD\n")

# Adjust prognostic covariates
extended_covariates <- cbind(propensity, covariates)

# Define X_test
X_test <- covariates

# Transform the outcome
time <- log(time)
time <- time - mean(time)

# Fit a Causal Horseshoe Forest
fit <- CausalShrinkageForest(
  y = time,
  status = status,
  X_train_control = extended_covariates,
  X_train_treat = covariates,
  treatment_indicator_train = treatment,
  outcome_type = "right-censored",
  prior_type_control = "horseshoe",
  prior_type_treat = "horseshoe",
  local_hp_control = 0.05 / sqrt(200),
  local_hp_treat = 0.05 / sqrt(200),
  global_hp_control = 0.05 / sqrt(200),
  global_hp_treat = 0.05 / sqrt(200),
  store_posterior_sample = TRUE,
  timescale = "log",
  number_of_trees_treat = 200,
  N_post = 5000,
  N_burn = 5000
)

# Evaluate C-index if survival package is available
if (requireNamespace("survival", quietly = TRUE)) {
  predicted_survtime <- fit$train_predictions
  cindex_result <- survival::concordance(Surv(time = time, event = status) ~ predicted_survtime)
  c_index <- cindex_result$concordance
  cat("C-index:", c_index, "\n")
} else {
  cat("Package 'survival' not available. Skipping C-index computation.\n")
}

# Compute posterior ATE samples
ate_samples <- rowMeans(fit$train_predictions_sample_treat)
mean_ate <- mean(ate_samples)
ci_95 <- quantile(ate_samples, probs = c(0.025, 0.975))

cat("Posterior mean ATE:", round(mean_ate, 3), "\n")
cat("95% credible interval: [", round(ci_95[1], 3), ", ", round(ci_95[2], 3), "]\n", sep = "")

# Plot histogram of ATE
hist(ate_samples,
     breaks = 30,
     col = col_treated,
     freq = FALSE,
     border = "white",
     xlab = "Average Treatment Effect",
     ylab = NULL,
     main = "Posterior distribution of ATE"
)
abline(v = mean_ate, col = "darkorchid4", lwd = 2)
abline(v = ci_95, col = "darkorchid4", lty = 2, lwd = 2)
legend("topright",
       legend = c("Mean", "95% CI"),
       col = c("darkorchid4", "darkorchid4"),
       lty = c(1, 2),
       lwd = 2)

# Prepare CATE summary
posterior_matrix <- fit$train_predictions_sample_treat
posterior_mean <- colMeans(posterior_matrix)
posterior_ci <- apply(posterior_matrix, 2, quantile, probs = c(0.025, 0.975))

df_cate <- data.frame(
  mean = posterior_mean,
  lower = posterior_ci[1, ],
  upper = posterior_ci[2, ]
)

# Sort by mean CATE
df_cate_sorted <- df_cate[order(df_cate$mean), ]
n <- nrow(df_cate_sorted)

# Plot CATE estimates
plot(
  x = df_cate_sorted$mean,
  y = 1:n,
  type = "n",
  xlab = "CATE per patient (95% credible interval)",
  ylab = "Patient index (sorted)",
  main = "Posterior CATE estimates",
  xlim = range(df_cate_sorted$lower, df_cate_sorted$upper)
)
segments(
  x0 = df_cate_sorted$lower,
  x1 = df_cate_sorted$upper,
  y0 = 1:n,
  y1 = 1:n,
  col = col_treated
)
lines(df_cate_sorted$mean, 1:n, col = "darkorchid4", lwd = 2)
abline(v = 0, col = "black", lwd = 2)

# Plot sigma trace
plot(fit$sigma,
     type = "l",
     xlab = "Iteration",
     ylab = expression(sigma),
     main = "Traceplot of Sigma",
     col = rgb(0, 0.5, 0, 0.5)
)
