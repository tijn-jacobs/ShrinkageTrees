library(ShrinkageTrees)
library(grf)
library(stochtree)
library(SoftBart)
library(doParallel)
library(foreach)

data_gen_nonlinear <- function(
    n = 200,
    p = 100,
    sigma = 1,
    sparse_tau = FALSE
) {
  
  # 1. Covariates
  X <- matrix(runif(n * p), nrow = n, ncol = p)
  
  # 2. Propensity score (nonlinear but low-dimensional signal)
  e_x <- pnorm(-0.5 + 0.4 * X[, 1] - 0.1 * X[, 3] + 0.3 * X[, 5])
  
  A <- rbinom(n, 1, e_x)
  
  # 3. Prognostic function (nonlinear Friedman-type)
  f_x <- 10 * sin(pi * X[, 1] * X[, 2]) +
    20 * (X[, 3] - 0.5)^2 +
    10 * X[, 4] +
    5 * X[, 5]
  
  # 4. Treatment effect function
  tau_x <-  2 * sin(pi * X[, 1] * X[, 3]) +
    3 * (X[, 2] - 0.5)^2 +
    2 * X[, 4]
  
  # Optional: add sparse linear high-dimensional component
  if (sparse_tau) {
    s <- 0.05
    beta_tau <- rbinom(p, 1, s) * rnorm(p)
    tau_x <- tau_x + X %*% beta_tau
  }
  
  # 5. Outcome
  Y_true <- f_x + (A - 1/2) * tau_x
  sd_Y <- sd(Y_true)
  Y_true <- Y_true / sd_Y
  Y <- Y_true + rnorm(n, 0, sigma)
  
  list(
    X = X,
    A = A,
    Y = as.numeric(Y),
    tau_true = as.numeric(tau_x) / sd_Y,
    f_true = as.numeric(f_x) / sd_Y,
    e_true = as.numeric(e_x)
  )
}

compute_metrics_cf <- function(tau_hat, tau_se, tau_true) {
  
  pehe <- sqrt(mean((tau_hat - tau_true)^2))
  
  ci_lower <- tau_hat - 1.96 * tau_se
  ci_upper <- tau_hat + 1.96 * tau_se
  
  coverage <- mean(tau_true >= ci_lower & tau_true <= ci_upper)
  length <- mean(ci_upper - ci_lower)
  
  c(
    PEHE = pehe,
    Coverage = coverage,
    Length = length
  )
}

fit_cf <- function(data, num.trees = 4000) {
  
  cf <- grf::causal_forest(
    X = data$X,
    Y = data$Y,
    W = data$A,
    num.trees = num.trees
  )
  
  pred <- predict(cf, estimate.variance = TRUE)
  
  list(
    tau_hat = as.numeric(pred$predictions),
    tau_se  = sqrt(as.numeric(pred$variance.estimates)),
    ate_hat = mean(pred$predictions)
  )
}

run_one_sim <- function(
    seed,
    n = 100,
    p = 1000,
    sigma = 1,
    num.trees = 4000
) {
  
  set.seed(seed)
  
  dat <- data_gen_nonlinear(
    n = n,
    p = p,
    sigma = sigma
  )
  
  cf_res <- fit_cf(dat, num.trees = num.trees)
  
  cf_metrics <- compute_metrics_cf(
    tau_hat = cf_res$tau_hat,
    tau_se  = cf_res$tau_se,
    tau_true = dat$tau_true
  )
  
  # Return a 1x3 matrix with rowname "CF"
  rbind(CF = cf_metrics)
}

run_simulation_parallel <- function(
    n_rep = 50,
    n = 100,
    p = 1000,
    sigma = 1,
    num.trees = 4000
) {
  
  # Each worker returns a matrix; wrap in list() so .combine = c makes a flat list
  results <- foreach(
    i = 1:n_rep,
    .combine = c,
    .packages = c("grf"),
    .export = c(
      "run_one_sim",
      "data_gen_nonlinear",
      "compute_metrics_cf",
      "fit_cf"
    )
  ) %dopar% {
    list(
      run_one_sim(
        seed = i,
        n = n,
        p = p,
        sigma = sigma,
        num.trees = num.trees
      )
    )
  }
  
  results
}

average_over_replications <- function(sim_list) {
  
  # sim_list is a list of matrices (each methods x metrics)
  methods <- rownames(sim_list[[1]])
  metrics <- colnames(sim_list[[1]])
  
  summary_mat <- matrix(
    NA_real_,
    nrow = length(methods),
    ncol = length(metrics),
    dimnames = list(methods, metrics)
  )
  
  for (m in methods) {
    mat_m <- do.call(
      rbind,
      lapply(sim_list, function(x) x[m, , drop = FALSE])
    )
    summary_mat[m, ] <- colMeans(mat_m)
  }
  
  summary_mat
}

# Retrieve command-line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) > 0) {
  num_cores <- as.integer(args[1]) - 1
} else {
  num_cores <- parallel::detectCores() - 1
}

cl <- makeCluster(num_cores)
registerDoParallel(cl)

cat("Number of cores being used (1 free):", num_cores, "\n")
cat("SIMULATION: revision continuous \n")

M <- 1000
n_train <- 100
sigma <- sqrt(3)

# RUN THE RESULTS HERE: three times p = 500, 1000, 5000
p_values <- c(500, 1000, 5000)

combined_results <- list()

for (p_val in p_values) {
  
  cat("Running simulation for p =", p_val, "\n")
  
  sim_res <- run_simulation_parallel(
    n_rep  = M,
    n      = n_train,
    p      = p_val,
    sigma  = sigma,
    num.trees = 4000
  )
  
  combined_results[[paste0("p_", p_val)]] <- list(
    raw = sim_res,
    summary = average_over_replications(sim_res)
  )
}

stopCluster(cl)

for (p_name in names(combined_results)) {
  cat("\nAveraged metrics for", p_name, "\n")
  print(combined_results[[p_name]]$summary)
}

# Define output file path
# (! NAME MUST BE FILENAME_output.rds !)
output_file <- file.path(Sys.getenv("TMPDIR"), "revision_continuous_output.rds")
cat("Saving all settings results to:", output_file, "\n")
saveRDS(combined_results, file = output_file)
cat("All results successfully saved in one file.\n")
