# ============================================================
# ACIC 2016 — Main simulation driver
#   - Serial over settings (1..77)
#   - Parallel over replications (1..M)
#   - Uses helper file defining: methods, run_one_rep, etc.
# ============================================================

suppressPackageStartupMessages({
  library(aciccomp2016)
  library(foreach)
  library(doParallel)
  library(SoftBart)
})

# ------------------------------------------------------------
# Source helpers
# ------------------------------------------------------------
source("/home/tjacobs/acic_helpers.R") 

# ------------------------------------------------------------
# Configuration
# ------------------------------------------------------------
cfg <- list(
  ndpost = 5,
  nskip  = 5,
  alpha = 0.05,
  seed_base = 1L,
  M = 100
)

# ------------------------------------------------------------
# Parallel backend
# ------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
num_cores <- if (length(args) > 0) as.integer(args[1]) - 1L else
  parallel::detectCores() - 1L
num_cores <- max(1L, num_cores)

cat("ACIC 2016 simulation\n")
cat("Cores used:", num_cores, "\n")
cat("Methods:", paste(names(methods), collapse = ", "), "\n")
cat("Replications per setting:", cfg$M, "\n\n")

registerDoParallel(cores = 5)

# ------------------------------------------------------------
# Fixed covariates (numeric design matrix)
# ------------------------------------------------------------
X <- model.matrix(~ . - 1, data = input_2016)

# ------------------------------------------------------------
# Main loop: serial over settings, parallel over reps
# ------------------------------------------------------------
all_results <- vector("list", 77)

for (setting in 1:77) {
  cat(format(Sys.time()), "Starting setting", setting, "of 77\n")
  
  res_setting <- foreach(
    rep = 1:cfg$M,
    .combine  = rbind,
    .packages = c("aciccomp2016", "BART", "SoftBart"),
    .export   = c("run_one_rep", "methods", "cfg", "X")
  ) %dopar% {
    run_one_rep(setting = setting, rep = rep, X = X, cfg = cfg, methods = methods)
  }
  
  all_results[[setting]] <- res_setting
  cat(format(Sys.time()), "Finished setting", setting, "\n\n")
}

# ------------------------------------------------------------
# Combine results
# ------------------------------------------------------------
results <- do.call(rbind, all_results)

# ------------------------------------------------------------
# Save results
# ------------------------------------------------------------
output_file <- file.path(Sys.getenv("TMPDIR"), "acic_new_output.rds")

saveRDS(
  list(
    config = cfg,
    methods = names(methods),
    results = results,
    all_results = all_results
  ),
  file = output_file
)

cat("All results saved to:", output_file, "\n")























## =========================================================
## Load ACIC results and inspect / summarize them
## =========================================================

# Load results
saved <- readRDS(
  "/Users/tijnjacobs/Library/CloudStorage/OneDrive-VrijeUniversiteitAmsterdam/Documents/GitHub/ShrinkageTrees/revision/acic/acic_new_output.rds"
)

config  <- saved$config
results <- saved$results

summary_overall <- data.frame(
  n_datasets = nrow(results),
  mean_pehe  = mean(results$pehe),
  median_pehe = median(results$pehe),
  rmse_ate   = sqrt(mean(results$ate_error^2)),
  mean_ate_coverage = mean(results$ate_coverage),
  mean_ate_intlen   = mean(results$ate_intlen),
  mean_ite_coverage = mean(results$ite_coverage, na.rm = TRUE),
  mean_ite_intlen   = mean(results$ite_intlen, na.rm = TRUE)
)


summary_by_setting <- aggregate(
  cbind(pehe, ate_error, ite_coverage, ite_intlen) ~ setting,
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
  mean_cate_coverage = mean(results$ite_coverage, na.rm = TRUE),
  mean_cate_intlen   = mean(results$ite_intlen, na.rm = TRUE)
)

overall
overall_check

