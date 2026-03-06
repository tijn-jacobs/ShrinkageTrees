# ============================================================
# ACIC 2016 — Causal Forest (GRF)
#   - Serial over settings (1..77)
#   - Parallel over replications
#   - Uses GRF asymptotic inference (no posterior draws)
# ============================================================

suppressPackageStartupMessages({
  library(aciccomp2016)
  library(foreach)
  library(doParallel)
  library(grf)
})

# ------------------------------------------------------------
# Metric computation (GRF-style inference)
# ------------------------------------------------------------
compute_metrics_cf <- function(tau_hat, tau_se, tau_true, alpha = 0.05) {
  
  # PEHE
  pehe <- sqrt(mean((tau_hat - tau_true)^2))
  
  # ATE error
  ate_error <- mean(tau_hat) - mean(tau_true)
  
  # CATE confidence intervals
  z <- qnorm(1 - alpha / 2)
  lower <- tau_hat - z * tau_se
  upper <- tau_hat + z * tau_se
  
  cate_coverage <- mean(tau_true >= lower & tau_true <= upper)
  cate_intlen   <- mean(upper - lower)
  
  data.frame(
    pehe           = pehe,
    ate_error      = ate_error,
    cate_coverage  = cate_coverage,
    cate_intlen    = cate_intlen
  )
}

# ------------------------------------------------------------
# Fit causal forest
# ------------------------------------------------------------
fit_cf <- function(data, seed) {
  
  set.seed(seed)
  
  cf <- causal_forest(
    X = data$X,
    Y = data$Y,
    W = data$Z,
    W.hat = data$propensity,  # oracle PS (ACIC)
    num.trees = 4000,
    honesty = TRUE,
    seed = seed
  )
  
  pred <- predict(cf, estimate.variance = TRUE)
  
  list(
    tau_hat = as.numeric(pred$predictions),
    tau_se  = sqrt(pred$variance.estimates)
  )
}

# ------------------------------------------------------------
# One replication
# ------------------------------------------------------------
run_one_rep <- function(setting, rep, X, cfg) {
  
  sim <- dgp_2016(
    x = X,
    parameters  = setting,
    random.seed = rep
  )
  
  data <- list(
    X = X,
    Z = as.numeric(sim$z),
    Y = as.numeric(sim$y),
    propensity = as.numeric(sim$e),
    tau_true = as.numeric(sim$mu.1) - as.numeric(sim$mu.0)
  )
  
  seed <- cfg$seed_base + 1000L * setting + rep
  
  fit <- fit_cf(data, seed)
  
  metrics <- compute_metrics_cf(
    tau_hat  = fit$tau_hat,
    tau_se   = fit$tau_se,
    tau_true = data$tau_true,
    alpha    = cfg$alpha
  )
  
  cbind(
    setting = setting,
    rep     = rep,
    method  = "GRF",
    metrics
  )
}

# ------------------------------------------------------------
# Configuration
# ------------------------------------------------------------
cfg <- list(
  alpha = 0.05,
  seed_base = 1L,
  M = 100
)

# ------------------------------------------------------------
# Parallel backend
# ------------------------------------------------------------
n_cores <- 190
registerDoParallel(cores = n_cores)

cat("ACIC 2016 — GRF simulation\n")
cat("Cores used:", n_cores, "\n")
cat("Replications per setting:", cfg$M, "\n\n")

# ------------------------------------------------------------
# Fixed covariates
# ------------------------------------------------------------
X <- model.matrix(~ . - 1, data = input_2016)

# ------------------------------------------------------------
# Main loop
# ------------------------------------------------------------
no_set <- 5
all_results <- vector("list", no_set)

for (setting in 1:no_set) {
  
  cat(format(Sys.time()), "Starting setting", setting, "of 77\n")
  
  res_setting <- foreach(
    rep = 1:cfg$M,
    .combine  = rbind,
    .packages = "grf"
  ) %dopar% {
    run_one_rep(setting, rep, X, cfg)
  }
  
  all_results[[setting]] <- res_setting
  
  cat(format(Sys.time()), "Finished setting", setting, "\n\n")
}

# ------------------------------------------------------------
# Combine and save
# ------------------------------------------------------------
results <- do.call(rbind, all_results)

output_file <- file.path(Sys.getenv("TMPDIR"), "acic_grf_output.rds")

saveRDS(
  list(
    config  = cfg,
    method  = "GRF",
    results = results
  ),
  file = output_file
)

cat("All results saved to:", output_file, "\n")





















## =========================================================
## Load ACIC results and inspect / summarize them
## =========================================================

# Load results
saved <- readRDS(
  "/Users/tijnjacobs/Library/CloudStorage/OneDrive-VrijeUniversiteitAmsterdam/Documents/GitHub/ShrinkageTrees/revision/acic/acic_grf_output.rds"
)

# Saved_y <- saved
# saved_mu5 <- saved
saved_mu77 <- saved


config  <- saved$config
results <- saved$results

summary_overall <- data.frame(
  n_datasets = nrow(results),
  mean_pehe  = mean(results$pehe),
  median_pehe = median(results$pehe),
  rmse_ate   = sqrt(mean(results$ate_error^2)),
  mean_cate_coverage = mean(results$cate_coverage, na.rm = TRUE),
  mean_cate_intlen   = mean(results$cate_intlen, na.rm = TRUE)
)


summary_by_setting <- aggregate(
  cbind(pehe, ate_error, cate_coverage, cate_intlen) ~ setting,
  data = results,
  FUN = function(x) c(mean = mean(x), sd = sd(x))
)


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
  mean_cate_coverage = mean(results$cate_coverage, na.rm = TRUE),
  mean_cate_intlen   = mean(results$cate_intlen, na.rm = TRUE)
)

overall
overall_check

