# ============================================================
# ACIC simulation — Helper functions
# ============================================================

suppressPackageStartupMessages({
  library(BART)
})

# ------------------------------------------------------------
# Generic metric computation
# ------------------------------------------------------------
compute_metrics <- function(tau_hat, tau_draws, tau_true, alpha = 0.05) {
  
  # PEHE
  pehe <- sqrt(mean((tau_hat - tau_true)^2))
  
  # ATE error
  ate_error <- mean(tau_hat) - mean(tau_true)
  
  # ATE interval
  ate_draws <- rowMeans(tau_draws)
  ate_ci <- quantile(ate_draws, probs = c(alpha / 2, 1 - alpha / 2))
  ate_coverage <- mean(tau_true) >= ate_ci[1] &&
    mean(tau_true) <= ate_ci[2]
  ate_intlen <- unname(diff(ate_ci))
  
  # CATE intervals
  lower <- apply(tau_draws, 2, quantile, probs = alpha / 2)
  upper <- apply(tau_draws, 2, quantile, probs = 1 - alpha / 2)
  
  cate_coverage <- mean(tau_true >= lower & tau_true <= upper)
  cate_intlen   <- mean(upper - lower)
  
  data.frame(
    pehe          = pehe,
    ate_error     = ate_error,
    ate_coverage  = ate_coverage,
    ate_intlen    = ate_intlen,
    cate_coverage  = cate_coverage,
    cate_intlen    = cate_intlen
  )
}

# ------------------------------------------------------------
# Methods
# ------------------------------------------------------------
fit_bart <- function(data, ndpost, nskip, seed) {
  
  set.seed(seed)
  
  X <- data$X
  Z <- data$Z
  Y <- data$Y
  propensity <- data$propensity
  
  X_train <- cbind(X, propensity, Z = Z)
  
  fit <- wbart(
    x.train = X_train,
    y.train = Y,
    ndpost  = ndpost,
    nskip   = nskip,
    printevery = 0
  )
  
  X1 <- cbind(X, propensity, Z = 1)
  X0 <- cbind(X, propensity, Z = 0)
  
  mu1 <- predict(fit, X1)
  mu0 <- predict(fit, X0)
  
  tau_draws <- mu1 - mu0
  tau_hat   <- colMeans(tau_draws)
  
  list(tau_hat = tau_hat, tau_draws = tau_draws)
}

fit_dart <- function(data, ndpost, nskip, seed) {
  
  set.seed(seed)
  
  X <- data$X
  Z <- data$Z
  Y <- data$Y
  propensity <- data$propensity
  
  X_train <- cbind(X, propensity, Z = Z)
  
  fit <- wbart(
    x.train = X_train,
    y.train = Y,
    ndpost  = ndpost,
    nskip   = nskip,
    sparse  = TRUE,
    printevery = 0
  )
  
  X1 <- cbind(X, propensity, Z = 1)
  X0 <- cbind(X, propensity, Z = 0)
  
  mu1 <- predict(fit, X1)
  mu0 <- predict(fit, X0)
  
  tau_draws <- mu1 - mu0
  tau_hat   <- colMeans(tau_draws)
  
  list(tau_hat = tau_hat, tau_draws = tau_draws)
}

fit_softbart <- function(data, ndpost, nskip, seed) {
  
  set.seed(seed)
  
  X <- data$X
  Z <- data$Z
  Y <- data$Y
  propensity <- data$propensity
  
  # Training design
  X_train <- cbind(X, propensity, Z = Z)
  
  # Counterfactual test designs
  X1 <- cbind(X, propensity, Z = 1)
  X0 <- cbind(X, propensity, Z = 0)
  
  # SoftBART options
  opts <- SoftBart::Opts(
    num_burn = nskip,
    num_save = ndpost
  )
  
  hypers <- SoftBart::Hypers(
    X = X_train,
    Y = Y,
    num_tree = 50,
    temperature = 1
  )
  
  fit <- SoftBart::softbart(
    X      = X_train,
    Y      = Y,
    X_test = rbind(X1, X0),
    hypers = hypers,
    opts   = opts,
    verbose = FALSE
  )
  
  # Extract posterior draws
  n <- nrow(X)
  mu1 <- fit$y_hat_test[, 1:n]
  mu0 <- fit$y_hat_test[, (n + 1):(2 * n)]
  
  tau_draws <- mu1 - mu0
  tau_hat   <- colMeans(tau_draws)
  
  list(
    tau_hat   = tau_hat,
    tau_draws = tau_draws
  )
}



methods <- list(
  BART = function(data, cfg, seed) {
    fit_bart(
      data = data,
      ndpost = cfg$ndpost,
      nskip  = cfg$nskip,
      seed   = seed
    )
  },
  DART = function(data, cfg, seed) {
    fit_dart(
      data = data,
      ndpost = cfg$ndpost,
      nskip  = cfg$nskip,
      seed   = seed
    )
  },
  SoftBART = function(data, cfg, seed) {
    fit_softbart(
      data   = data,
      ndpost = cfg$ndpost,
      nskip  = cfg$nskip,
      seed   = seed
    )
  }
)


run_one_rep <- function(setting, rep, X, cfg, methods) {
  
  sim <- dgp_2016(
    x = X,
    parameters  = setting,
    random.seed = rep
  )
  
  sim_data <- list(
    X = X,
    Z = as.numeric(sim$z),
    Y = as.numeric(sim$y),
    propensity = as.numeric(sim$e),
    tau_true = as.numeric(sim$y.1) - as.numeric(sim$y.0)
  )
  
  out <- lapply(names(methods), function(m) {
    
    seed <- cfg$seed_base + 1000L * setting + rep
    fit  <- methods[[m]](sim_data, cfg, seed)
    
    metrics <- compute_metrics(
      tau_hat   = fit$tau_hat,
      tau_draws = fit$tau_draws,
      tau_true  = sim_data$tau_true,
      alpha     = cfg$alpha
    )
    
    cbind(
      setting = setting,
      rep     = rep,
      method  = m,
      metrics
    )
  })
  
  do.call(rbind, out)
}
