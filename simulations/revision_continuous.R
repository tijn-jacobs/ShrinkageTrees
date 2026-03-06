library(ShrinkageTrees)
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
  Y_true <- Y_true/sd_Y
  Y <- Y_true + rnorm(n, 0, sigma)
  
  list(
    X = X,
    A = A,
    Y = as.numeric(Y),
    tau_true = as.numeric(tau_x)/sd_Y,
    f_true = as.numeric(f_x)/sd_Y,
    e_true = as.numeric(e_x)
  )
}

compute_metrics <- function(tau_hat, tau_post, tau_true) {
  
  pehe <- sqrt(mean((tau_hat - tau_true)^2))
  
  ci_lower <- apply(tau_post, 2, quantile, 0.025)
  ci_upper <- apply(tau_post, 2, quantile, 0.975)
  
  coverage <- mean(tau_true >= ci_lower &
                     tau_true <= ci_upper)
  
  length <- mean(ci_upper - ci_lower)
  
  c(PEHE = pehe,
    Coverage = coverage,
    Length = length)
}

compute_metrics_cf <- function(tau_hat, tau_se, tau_true) {
  
  pehe <- sqrt(mean((tau_hat - tau_true)^2))
  
  ci_lower <- tau_hat - 1.96 * tau_se
  ci_upper <- tau_hat + 1.96 * tau_se
  
  coverage <- mean(tau_true >= ci_lower &
                     tau_true <= ci_upper)
  
  length <- mean(ci_upper - ci_lower)
  
  c(PEHE = pehe,
    Coverage = coverage,
    Length = length)
}



fit_chf <- function(
    data,
    k = 0.1,
    N_post = 2500,
    N_burn = 2500,
    m = 200
) {
  
  e_hat <- as.numeric(data$e_true)
  sd_Y <- 1#sd(data$Y)
  Y_train <- data$Y/sd_Y
  
  chf <- ShrinkageTrees::CausalShrinkageForest(
    y = Y_train,
    X_train_control = cbind(data$X, e_hat),
    X_train_treat   = data$X,
    treatment_indicator_train = data$A,
    outcome_type = "continuous",

    number_of_trees_control = m,
    number_of_trees_treat   = m,
    
    prior_type_control = "horseshoe",
    prior_type_treat   = "horseshoe",
    
    local_hp_control  = k / sqrt(m),
    local_hp_treat    = k / sqrt(m),
    global_hp_control = k / sqrt(m),
    global_hp_treat   = k / sqrt(m),
    
    N_post = N_post,
    N_burn = N_burn,
    store_posterior_sample = TRUE,
    verbose = FALSE
  )
  
  list(
    tau_hat  = as.numeric(chf$train_predictions_treat)*sd_Y,
    tau_post = chf$train_predictions_sample_treat*sd_Y,
    ate_hat  = mean(chf$train_predictions_treat)*sd_Y,
    ate_post = rowMeans(chf$train_predictions_sample_treat)*sd_Y
  )
}

fit_bart_s <- function(
    data,
    N_post = 2500,
    N_burn = 2500
) {
  
  X_train <- cbind(data$X, A = data$A, e = data$e_true)
  
  # Counterfactual test matrix
  X1 <- cbind(data$X, A = 1, e = data$e_true)
  X0 <- cbind(data$X, A = 0, e = data$e_true)
  X_test_full <- rbind(X1, X0)
  
  bart_fit <- stochtree::bart(
    X_train = X_train,
    y_train = data$Y,
    X_test  = X_test_full,
    num_gfr = 0,
    num_burnin = N_burn,
    num_mcmc = N_post,
    general_params = list(verbose = FALSE)
  )
  
  preds <- bart_fit$y_hat_test  # [draws × 2n]
  
  n <- nrow(data$X)
  
  pred1 <- preds[, 1:n]
  pred0 <- preds[, (n + 1):(2 * n)]
  
  tau_post <- pred1 - pred0
  tau_hat  <- colMeans(tau_post)
  
  list(
    tau_hat  = as.numeric(tau_hat),
    tau_post = tau_post,
    ate_hat  = mean(tau_hat),
    ate_post = rowMeans(tau_post)
  )
}

fit_bart <- function(
    data,
    N_post = 2500,
    N_burn = 2500
) {
  
  X_train <- cbind(data$X, A = data$A, e = data$e_true)
  
  X1 <- cbind(data$X, A = 1, e = data$e_true)
  X0 <- cbind(data$X, A = 0, e = data$e_true)
  X_test_full <- rbind(X1, X0)
  
  bart_fit <- capture.output(
    fit <- BART::wbart(
      x.train = X_train,
      y.train = data$Y,
      x.test  = X_test_full,
      nskip   = N_burn,
      ndpost  = N_post,
      printevery = 999999999
    ),
    type = "output"
  )
  bart_fit <- fit
  
  preds <- bart_fit$yhat.test  # [draws × 2n]
  
  n <- nrow(data$X)
  
  pred1 <- preds[, 1:n]
  pred0 <- preds[, (n + 1):(2 * n)]
  
  tau_post <- pred1 - pred0
  tau_hat  <- colMeans(tau_post)
  
  list(
    tau_hat  = as.numeric(tau_hat),
    tau_post = tau_post,
    ate_hat  = mean(tau_hat),
    ate_post = rowMeans(tau_post)
  )
}

fit_dart <- function(
    data,
    N_post = 2500,
    N_burn = 2500
) {
  
  X_train <- cbind(data$X, A = data$A, e = data$e_true)
  
  X1 <- cbind(data$X, A = 1, e = data$e_true)
  X0 <- cbind(data$X, A = 0, e = data$e_true)
  X_test_full <- rbind(X1, X0)
  
  dart_fit <- capture.output(
    fit <- BART::wbart(
      x.train = X_train,
      y.train = data$Y,
      x.test  = X_test_full,
      sparse  = TRUE,
      nskip   = N_burn,
      ndpost  = N_post,
      printevery = 999999999
    ),
    type = "output"
  )
  dart_fit <- fit
  
  preds <- dart_fit$yhat.test  # [draws × 2n]
  
  n <- nrow(data$X)
  
  pred1 <- preds[, 1:n]
  pred0 <- preds[, (n + 1):(2 * n)]
  
  tau_post <- pred1 - pred0
  tau_hat  <- colMeans(tau_post)
  
  list(
    tau_hat  = as.numeric(tau_hat),
    tau_post = tau_post,
    ate_hat  = mean(tau_hat),
    ate_post = rowMeans(tau_post)
  )
}



fit_bcf <- function(
    data,
    m_mu = 200,
    m_tau = 50,
    N_post = 2500,
    N_burn = 2500
) {
  
  bcf_fit <- stochtree::bcf(
    X_train = data$X,
    y_train = data$Y,
    Z_train = data$A,
    propensity_train = data$e_true,
    num_gfr = 0,
    num_burnin = N_burn,
    num_mcmc = N_post,
    general_params = list(
      verbose = FALSE
    )
  )
  
  # Posterior tau draws
  # stochtree stores tau draws inside the model object
  tau_post <- bcf_fit$tau_hat
  
  # Dimensions: [draws × n]
  tau_hat <- colMeans(tau_post)
  
  list(
    tau_hat  = as.numeric(tau_hat),
    tau_post = tau_post,
    ate_hat  = mean(tau_hat),
    ate_post = rowMeans(tau_post)
  )
}

fit_softbart <- function(
    data,
    m = 200,
    N_post = 2500,
    N_burn = 2500
) {
  
  # Design matrix (S-learner)
  X_train <- cbind(data$X, A = data$A, e = data$e_true)
  
  # Counterfactual matrices
  X1 <- cbind(data$X, A = 1, e = data$e_true)
  X0 <- cbind(data$X, A = 0, e = data$e_true)
  
  # Combine into one test matrix to avoid refitting twice
  X_test_full <- rbind(X1, X0)
  
  # Hyperparameters
  hypers <- SoftBart::Hypers(
    X = X_train,
    Y = data$Y,
    num_tree = m
  )
  
  opts <- SoftBart::Opts(
    num_burn = N_burn,
    num_save = N_post
  )
  
  sb_fit <- SoftBart::softbart(
    X = X_train,
    Y = data$Y,
    X_test = X_test_full,
    hypers = hypers,
    opts = opts,
    verbose = FALSE
  )
  
  # Extract posterior predictions
  # y_hat_test: [num_save × (2n)]
  preds <- sb_fit$y_hat_test
  
  n <- nrow(data$X)
  
  pred1 <- preds[, 1:n]
  pred0 <- preds[, (n + 1):(2 * n)]
  
  tau_post <- pred1 - pred0
  tau_hat  <- colMeans(tau_post)
  
  list(
    tau_hat  = as.numeric(tau_hat),
    tau_post = tau_post,
    ate_hat  = mean(tau_hat),
    ate_post = rowMeans(tau_post)
  )
}


run_one_sim <- function(
    seed,
    n = 100,
    p = 1000,
    sigma = 1,
    N_post = 2000,
    N_burn = 2000
) {
  
  set.seed(seed)
  
  dat <- data_gen_nonlinear(
    n = n,
    p = p,
    sigma = sigma
  )
  
  # Fit models
  chf0.1_res   <- fit_chf(dat, k = 0.1, N_post = N_post, N_burn = N_burn)
  chf0.2_res   <- fit_chf(dat, k = 0.2, N_post = N_post, N_burn = N_burn)
  chf0.3_res   <- fit_chf(dat, k = 0.3, N_post = N_post, N_burn = N_burn)
  #bart_s_res   <- fit_bart_s(dat, N_post = N_post, N_burn = N_burn)
  bart_res     <- fit_bart(dat, N_post = N_post, N_burn = N_burn)
  dart_res     <- fit_dart(dat, N_post = N_post, N_burn = N_burn)
  bcf_res      <- fit_bcf(dat, N_post = N_post, N_burn = N_burn)
  softbart_res <- fit_softbart(dat, N_post = N_post, N_burn = N_burn)
  
  # Metrics
  chf0.1_metrics      <- compute_metrics(chf0.1_res$tau_hat,
                                      chf0.1_res$tau_post,
                                      dat$tau_true)
  
  chf0.2_metrics      <- compute_metrics(chf0.2_res$tau_hat,
                                      chf0.2_res$tau_post,
                                      dat$tau_true)
  
  chf0.3_metrics      <- compute_metrics(chf0.3_res$tau_hat,
                                      chf0.3_res$tau_post,
                                      dat$tau_true)
  
  # bart_s_metrics   <- compute_metrics(bart_s_res$tau_hat,
  #                                     bart_s_res$tau_post,
  #                                     dat$tau_true)
  
  bart_metrics     <- compute_metrics(bart_res$tau_hat,
                                      bart_res$tau_post,
                                      dat$tau_true)
  
  dart_metrics     <- compute_metrics(dart_res$tau_hat,
                                      dart_res$tau_post,
                                      dat$tau_true)
  
  bcf_metrics      <- compute_metrics(bcf_res$tau_hat,
                                      bcf_res$tau_post,
                                      dat$tau_true)
  
  softbart_metrics <- compute_metrics(softbart_res$tau_hat,
                                      softbart_res$tau_post,
                                      dat$tau_true)
  
  rbind(
    CHF0.1       = chf0.1_metrics,
    CHF0.2       = chf0.2_metrics,
    CHF0.3       = chf0.3_metrics,
    #BART_s    = bart_s_metrics,
    BART      = bart_metrics,
    DART      = dart_metrics,
    BCF       = bcf_metrics,
    SoftBART  = softbart_metrics
  )
}


run_simulation_parallel <- function(
    n_rep = 50,
    n = 100,
    p = 1000,
    sigma = 1,
    N_post = 2000,
    N_burn = 2000,
    n_cores = parallel::detectCores() - 1
) {
  
  results <- foreach(
    i = 1:n_rep,
    .combine = "list",
    .packages = c("ShrinkageTrees",
                  "BART",
                  "stochtree",
                  "SoftBart"),
    .export = c(
      "run_one_sim",
      "data_gen_nonlinear",
      "compute_metrics",
      "fit_chf",
      "fit_bart",
      "fit_dart",
      "fit_bcf",
      "fit_softbart"
    )
  ) %dopar% {
    
    run_one_sim(
      seed = i,
      n = n,
      p = p,
      sigma = sigma,
      N_post = N_post,
      N_burn = N_burn
    )
  }
  
  results
}






# Retrieve command-line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) > 0) {
  num_cores <- as.integer(args[1]) - 1
} else {
  num_cores <- parallel::detectCores() - 1
}

registerDoParallel(cores = num_cores)

cat("Number of cores being used (1 free):", num_cores, "\n")
cat("SIMULATION: revision continuous \n")

M <- num_cores
n_train <- 100
sigma <- sqrt(3)
N_post <- 50
N_burn <- 25

# RUN THE RESULTS HERE: three times p = 500, 1000, 50000
p_values <- c(500, 1000, 5000)  # adjust if needed

combined_results <- list()

for (p_val in p_values) {
  
  cat("Running simulation for p =", p_val, "\n")
  
  sim_res <- run_simulation_parallel(
    n_rep  = M,
    n      = n_train,
    p      = p_val,
    sigma  = sigma,
    N_post = N_post,
    N_burn = N_burn,
    n_cores = num_cores
  )
  
  # Aggregate across replications
  methods <- rownames(sim_res[[1]])
  
  summary_mat <- matrix(NA,
                        nrow = length(methods),
                        ncol = 3,
                        dimnames = list(methods,
                                        c("PEHE", "Coverage", "Length")))
  
  for (m in methods) {
    mat <- do.call(rbind,
                   lapply(sim_res, function(x) x[m, ]))
    summary_mat[m, ] <- colMeans(mat)
  }
  
  combined_results[[paste0("p_", p_val)]] <- list(
    raw = sim_res,
    summary = summary_mat
  )
}

# Define output file path 
# (! NAME MUST BE FILENAME_output.rds !)
output_file <- file.path(Sys.getenv("TMPDIR"), "revision_continuous_output.rds")
cat("Saving all settings results to:", output_file, "\n")
saveRDS(combined_results, file = output_file)
cat("All results successfully saved in one file.\n")

























































# 
# N_post <- 100
# N_burn <- 100
# #dat <- data_gen_nonlinear(n = 100, p = 1000, sigma = 1)
# 
# # Fit the models
# cf_res  <- fit_cf(dat)
# chf_res <- fit_chf(dat, k = 0.2, N_post = N_post, N_burn = N_burn)
# bart_s_res <- fit_bart_s(dat, N_post = N_post, N_burn = N_burn)
# bart_res <- fit_bart(dat, N_post = N_post, N_burn = N_burn)
# dart_res <- fit_dart(dat, N_post = N_post, N_burn = N_burn)
# bcf_res <- fit_bcf(dat, N_post = N_post, N_burn = N_burn)
# softbart_res <- fit_softbart(dat, N_post = N_post, N_burn = N_burn)
# 
# # CF metrics
# pehe_cf <- sqrt(mean((cf_res$tau_hat - dat$tau_true)^2))
# ci_lower_cf <- cf_res$tau_hat - 1.96 * cf_res$tau_se
# ci_upper_cf <- cf_res$tau_hat + 1.96 * cf_res$tau_se
# coverage_cf <- mean(dat$tau_true >= ci_lower_cf &
#                       dat$tau_true <= ci_upper_cf)
# length_cf <- mean(ci_upper_cf - ci_lower_cf)
# 
# # CHF metrics
# pehe_chf <- sqrt(mean((chf_res$tau_hat - dat$tau_true)^2))
# ci_lower_chf <- apply(chf_res$tau_post, 2, quantile, 0.025)
# ci_upper_chf <- apply(chf_res$tau_post, 2, quantile, 0.975)
# coverage_chf <- mean(dat$tau_true >= ci_lower_chf &
#                        dat$tau_true <= ci_upper_chf)
# length_chf <- mean(ci_upper_chf - ci_lower_chf)
# 
# # BART_s metrics
# pehe_bart_s <- sqrt(mean((bart_s_res$tau_hat - dat$tau_true)^2))
# ci_lower_bart_s <- apply(bart_s_res$tau_post, 2, quantile, 0.025)
# ci_upper_bart_s <- apply(bart_s_res$tau_post, 2, quantile, 0.975)
# coverage_bart_s <- mean(dat$tau_true >= ci_lower_bart_s &
#                         dat$tau_true <= ci_upper_bart_s)
# length_bart_s <- mean(ci_upper_bart_s - ci_lower_bart_s)
# 
# # BART metrics
# pehe_bart <- sqrt(mean((bart_res$tau_hat - dat$tau_true)^2))
# ci_lower_bart <- apply(bart_res$tau_post, 2, quantile, 0.025)
# ci_upper_bart <- apply(bart_res$tau_post, 2, quantile, 0.975)
# coverage_bart <- mean(dat$tau_true >= ci_lower_bart &
#                         dat$tau_true <= ci_upper_bart)
# length_bart <- mean(ci_upper_bart - ci_lower_bart)
# 
# # DART metrics
# pehe_dart <- sqrt(mean((dart_res$tau_hat - dat$tau_true)^2))
# ci_lower_dart <- apply(dart_res$tau_post, 2, quantile, 0.025)
# ci_upper_dart <- apply(dart_res$tau_post, 2, quantile, 0.975)
# coverage_dart <- mean(dat$tau_true >= ci_lower_dart &
#                         dat$tau_true <= ci_upper_dart)
# length_dart <- mean(ci_upper_dart - ci_lower_dart)
# 
# 
# # bcf metrics
# pehe_bcf <- sqrt(mean((bcf_res$tau_hat - dat$tau_true)^2))
# ci_lower_bcf <- apply(bcf_res$tau_post, 2, quantile, 0.025)
# ci_upper_bcf <- apply(bcf_res$tau_post, 2, quantile, 0.975)
# coverage_bcf <- mean(dat$tau_true >= ci_lower_bcf &
#                        dat$tau_true <= ci_upper_bcf)
# length_bcf <- mean(ci_upper_bcf - ci_lower_bcf)
# 
# # softbart metrics
# pehe_softbart <- sqrt(mean((softbart_res$tau_hat - dat$tau_true)^2))
# ci_lower_softbart <- apply(softbart_res$tau_post, 2, quantile, 0.025)
# ci_upper_softbart <- apply(softbart_res$tau_post, 2, quantile, 0.975)
# coverage_softbart <- mean(dat$tau_true >= ci_lower_softbart &
#                             dat$tau_true <= ci_upper_softbart)
# length_softbart <- mean(ci_upper_softbart - ci_lower_softbart)
# 
# 
# list(
#   CF = c(PEHE = pehe_cf,
#          Coverage = coverage_cf,
#          Length = length_cf),
#   
#   CHF = c(PEHE = pehe_chf,
#           Coverage = coverage_chf,
#           Length = length_chf),
#   
#   BART_s = c(PEHE = pehe_bart_s,
#           Coverage = coverage_bart_s,
#           Length = length_bart_s),
#   
#   BART = c(PEHE = pehe_bart,
#            Coverage = coverage_bart,
#            Length = length_bart),
#   
#   SoftBART = c(PEHE = pehe_softbart,
#            Coverage = coverage_softbart,
#            Length = length_softbart),
#   
#   DART = c(PEHE = pehe_dart,
#            Coverage = coverage_dart,
#            Length = length_dart),
#   
#   BCF = c(PEHE = pehe_bcf,
#            Coverage = coverage_bcf,
#            Length = length_bcf)
# )
# 
# 
# 




