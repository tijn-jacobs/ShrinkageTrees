# ============================================================
# ACIC 2016 — BART baseline using BART::wbart()
# Design:
#   - Serial loop over 77 settings
#   - Parallel loop over replications (M = 100)
#   - Robust progress logging
# ============================================================

suppressPackageStartupMessages({
  library(aciccomp2016)
  library(BART)
  library(foreach)
  library(doParallel)
})

# ------------------------------------------------------------
# Configuration
# ------------------------------------------------------------
cfg <- list(
  ndpost = 5000,
  nskip  = 5000,
  compute_intervals = TRUE,
  alpha = 0.05,
  seed_base = 1L,
  M = 100
)

# ------------------------------------------------------------
# Parallel backend
# ------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  n_cores <- as.integer(args[1]) - 1
} else {
  n_cores <- parallel::detectCores() - 1
}
n_cores <- max(1L, n_cores)

cl <- makeCluster(n_cores)
registerDoParallel(cl)

cat("ACIC 2016 BART simulation\n")
cat("Cores used:", n_cores, "\n")
cat("Replications per setting:", cfg$M, "\n\n")

# ------------------------------------------------------------
# Fixed covariates (numeric design matrix)
# ------------------------------------------------------------
X <- model.matrix(~ . - 1, data = input_2016)

# ------------------------------------------------------------
# Helper: Fit BART and compute ITEs
# ------------------------------------------------------------
fit_bart <- function(X, Z, Y, ndpost, nskip, seed) {
  
  set.seed(seed)
  
  X_train <- cbind(X, Z = Z)
  
  n <- nrow(X_train)
  p_noise <- 500
  X_noise <- matrix(runif(n*p_noise), nrow = n, ncol = p_noise)
  X_train <- cbind(X_train, X_noise)
  
  fit <- wbart(
    x.train = X_train,
    y.train = Y,
    ndpost  = ndpost,
    nskip   = nskip
  )
  
  # reuse the SAME noise for counterfactuals
  X1 <- cbind(X, Z = 1, X_noise)
  X0 <- cbind(X, Z = 0, X_noise)
  
  mu1 <- predict(fit, X1)
  mu0 <- predict(fit, X0)
  
  tau_draws <- mu1 - mu0
  tau_hat   <- colMeans(tau_draws)
  
  list(tau_hat = tau_hat, tau_draws = tau_draws)
}

# ------------------------------------------------------------
# Metrics
# ------------------------------------------------------------
metric_pehe <- function(tau_hat, tau_true) {
  sqrt(mean((tau_hat - tau_true)^2))
}

metric_ate_error <- function(tau_hat, tau_true) {
  mean(tau_hat) - mean(tau_true)
}

metric_cate_coverage <- function(tau_draws, tau_true, alpha) {
  lower <- apply(tau_draws, 2, quantile, probs = alpha / 2)
  upper <- apply(tau_draws, 2, quantile, probs = 1 - alpha / 2)
  mean(tau_true >= lower & tau_true <= upper)
}

metric_cate_intlen <- function(tau_draws, alpha) {
  lower <- apply(tau_draws, 2, quantile, probs = alpha / 2)
  upper <- apply(tau_draws, 2, quantile, probs = 1 - alpha / 2)
  mean(upper - lower)
}

metric_ate_interval <- function(tau_draws, alpha) {
  ate_draws <- rowMeans(tau_draws)
  quantile(ate_draws, probs = c(alpha / 2, 1 - alpha / 2))
}

metric_ate_coverage <- function(tau_draws, tau_true, alpha) {
  ate_true <- mean(tau_true)
  ci <- metric_ate_interval(tau_draws, alpha)
  ate_true >= ci[1] && ate_true <= ci[2]
}

metric_ate_intlen <- function(tau_draws, alpha) {
  ci <- metric_ate_interval(tau_draws, alpha)
  unname(ci[2] - ci[1])
}

# ------------------------------------------------------------
# One replication (unit of parallel work)
# ------------------------------------------------------------
run_one_rep <- function(setting, rep, X, cfg) {
  
  sim <- dgp_2016(
    x = X,
    parameters  = setting,
    random.seed = rep
  )
  
  Z <- as.numeric(sim$z)
  Y <- as.numeric(sim$y)
  tau_true <- as.numeric(sim$y.1) - as.numeric(sim$y.0)
  
  seed <- cfg$seed_base + 1000L * setting + rep
  
  fit <- fit_bart(
    X = X,
    Z = Z,
    Y = Y,
    ndpost = cfg$ndpost,
    nskip  = cfg$nskip,
    seed   = seed
  )
  
  pehe <- metric_pehe(fit$tau_hat, tau_true)
  ate_err <- metric_ate_error(fit$tau_hat, tau_true)
  
  ate_cov  <- metric_ate_coverage(fit$tau_draws, tau_true, cfg$alpha)
  ate_len  <- metric_ate_intlen(fit$tau_draws, cfg$alpha)
  
  if (cfg$compute_intervals) {
    ite_cov  <- metric_cate_coverage(fit$tau_draws, tau_true, cfg$alpha)
    ite_len <- metric_cate_intlen(fit$tau_draws, cfg$alpha)
  } else {
    ite_cov  <- NA_real_
    ite_len <- NA_real_
  }
  
  data.frame(
    setting = setting,
    rep     = rep,
    pehe    = pehe,
    ate_error = ate_err,
    ate_coverage = ate_cov,
    ate_intlen  = ate_len,
    ite_coverage = ite_cov,
    ite_intlen   = ite_len
  )
}
# ------------------------------------------------------------
# Main loop: serial over settings, parallel over reps
# ------------------------------------------------------------
all_results <- vector("list", 77)

for (setting in 1:77) {
  
  msg <- paste(format(Sys.time()),
               "Starting setting", setting, "of 77\n")
  cat(msg)
  
  res_setting <- foreach(
    rep = 1:cfg$M,
    .combine  = rbind,
    .packages = c("aciccomp2016", "BART")
  ) %dopar% {
    run_one_rep(setting, rep, X, cfg)
  }
  
  all_results[[setting]] <- res_setting
  
  msg <- paste(format(Sys.time()),
               "Finished setting", setting, "\n\n")
  cat(msg)
}

stopCluster(cl)

# ------------------------------------------------------------
# Combine and summarize
# ------------------------------------------------------------
results <- do.call(rbind, all_results)

summary_overall <- data.frame(
  n_datasets = nrow(results),
  mean_pehe  = mean(results$pehe),
  median_pehe = median(results$pehe),
  rmse_ate   = sqrt(mean(results$ate_error^2)),
  mean_ate_coverage = mean(results$ate_coverage),
  mean_ate_intlen   = mean(results$ate_intlen),
  mean_cate_coverage = mean(results$ite_coverage, na.rm = TRUE),
  mean_cate_intlen   = mean(results$ite_intlen, na.rm = TRUE)
)


summary_by_setting <- aggregate(
  cbind(pehe, ate_error, ite_coverage, ite_intlen) ~ setting,
  data = results,
  FUN = function(x) c(mean = mean(x), sd = sd(x))
)

# ------------------------------------------------------------
# Save results
# ------------------------------------------------------------
output_file <- file.path(Sys.getenv("TMPDIR"), "acic500_output.rds")

saveRDS(
  list(
    config = cfg,
    results = results,
    summary_overall = summary_overall,
    summary_by_setting = summary_by_setting
  ),
  file = output_file
)

cat("All results saved to:", output_file, "\n")









## =========================================================
## Load ACIC results and inspect / summarize them
## =========================================================

# Load results
saved <- readRDS(
  "/Users/tijnjacobs/Library/CloudStorage/OneDrive-VrijeUniversiteitAmsterdam/Documents/GitHub/ShrinkageTrees/simulations revision/acic500_output.rds"
)

# Unpack
config  <- saved$config
results <- saved$results
overall <- saved$summary_overall
by_set  <- saved$summary_by_setting

# ---------------------------------------------------------
# Basic inspection
# ---------------------------------------------------------
str(saved)
head(results)
tail(results)

# What was actually run?
config
length(unique(results$setting))
nrow(results)

# ---------------------------------------------------------
# Overall summaries (recomputed from raw results)
# ---------------------------------------------------------
overall_check <- data.frame(
  n_datasets = nrow(results),
  mean_pehe  = mean(results$pehe),
  median_pehe = median(results$pehe),
  rmse_ate   = sqrt(mean(results$ate_error^2)),
  mean_cate_coverage = mean(results$ite_coverage, na.rm = TRUE),
  mean_cate_intlen   = mean(results$ite_intlen, na.rm = TRUE)
)

overall
overall_check