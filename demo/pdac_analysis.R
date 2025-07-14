## Demo: Analysis of the TCGA PAAD dataset (pdac_plain)
## ShrinkageTrees package
##
## This demo loads the pdac dataset and shows the data analysis

# Load package
library(ShrinkageTrees)

# Load the data
data("pdac_plain")

# Retrieve the data
time <- pdac_plain$time
status <- pdac_plain$status
treatment <- pdac_plain$treatment
covariates <- pdac_plain[, !(colnames(pdac_plain) %in% 
                               c("time", "status", "treatment"))]

# Estimate the propensity scores
propensity_fit <- HorseTrees(y = treatment,
                             X_train = covariates,
                             outcome_type = "binary",
                             k = 0.1,
                             N_post = 1000,
                             N_burn = 1000
)

# Retrieve estimated propensity scores
propensity <- pnorm(propensity_fit$train_predictions)

# Ask user if they want to plot overlap
plot_overlap <- readline(prompt = "Plot propensity score overlap? (y/n): ")

if (tolower(plot_overlap) %in% c("y", "yes")) {
  
  # Define colors
  col_treated <- rgb(0, 0.5, 0, 0.5)  # deep green
  col_control <- rgb(1, 0.5, 0, 0.5)  # deep orange
  
  # Plot control group histogram
  hist(p0,
       breaks = 10,
       xlim = c(0.3, 0.7),
       col = col_control,
       xlab = "Propensity score",
       ylab = NULL,
       main = "Propensity score overlap"
  )
  
  # Add treated group histogram
  hist(p1,
       breaks = 10,
       xlim = c(0.3, 0.7),
       col = col_treated,
       add = TRUE
  )
  
  # Add legend
  legend("topright",
         legend = c("Control", "Treated"),
         fill = c(col_control, col_treated)
  )
  
} else {
  cat("Skipping overlap plot.\n")
}

# Placeholder for your analysis
cat("\n Analysis section TBD \n")

# Adjust the prognostic model also for the propensity score
extended_covariates <- cbind(propensity, covariates)

# Transform the outcome
time <- log(time)
time <- time - mean(time)
time <- time/sd(time)

# Fit a Causal Horshoe Forest
fit <- CausalHorseForest(
  y = time,
  status = status,
  X_control_train = extended_covariates,
  X_treat_train = covariates,
  X_control_test = extended_covariates,
  X_treat_test = covariates,
  treatment_indicator_train = treatment,
  treatment_indicator_test = rep(1, nrow(X_test)),
  store_posterior_sample_control = FALSE,
  store_posterior_sample_treat = TRUE,
  scale = "log",
  verbose = TRUE,
  alpha_local_control = 0.05/sqrt(200),
  alpha_global_control = 0.05/sqrt(200),
  alpha_local_treat = 0.05/sqrt(200),
  alpha_global_treat = 0.05/sqrt(200),
  omega_control = 1/2,
  omega_treat = 1/2,
  number_of_trees_treat = 200,
  N_post = 1000,
  N_burn = 1000
)

if (requireNamespace("survival", quietly = TRUE)) {
  # Compute C-index
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

# Create a data frame (optional, for further inspection or plotting)
df <- data.frame(ate = ate_samples)

# Present the results nicely
cat("Posterior mean ATE:", round(mean_ate, 3), "\n")
cat("95% credible interval: [", round(ci_95[1], 3), ", ", round(ci_95[2], 3), "]\n", sep = "")

# Plot histogram
hist(ate_samples,
     breaks = 30,
     col = col_treated,
     freq = FALSE,
     border = "white",
     xlab = "Average Treatment Effect",
     ylab = NULL,
     main = "Posterior distribution of ATE")
  
# Add vertical lines for mean and credible interval
abline(v = mean_ate, col = "darkorchid4", lwd = 2)
abline(v = ci_95, col = "darkorchid4", lty = 2, lwd = 2)
  
# Add legend
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


# Initialize empty plot
plot(
  x = df_cate_sorted$mean,
  y = 1:n,
  type = "n",
  xlab = "CATE per patient (95% credible interval)",
  ylab = "Patient index (sorted)",
  main = "Posterior CATE estimates",
  xlim = range(df_cate_sorted$lower, df_cate_sorted$upper)
)

# Draw intervals
segments(
  x0 = df_cate_sorted$lower,
  x1 = df_cate_sorted$upper,
  y0 = 1:n,
  y1 = 1:n,
  col = col_treated
)

# Draw posterior means
lines(df_cate_sorted$mean, 1:n, col = "darkorchid4", lwd = 2)

# Reference line at 0
abline(v = 0, col = "black", lwd = 2)


# Plot the sigma trace
plot(fit$sigma,
     type = "l",
     xlab = "Iteration",
     ylab = expression(sigma),
     main = "Traceplot of Sigma",
     col = rgb(0, 0.5, 0, 0.5)  
)

