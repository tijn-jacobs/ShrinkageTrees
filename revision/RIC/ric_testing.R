library(ShrinkageTrees)
library(AFTrees)
library(foreach)
library(doParallel)

source("/home/tjacobs/evaluation_functions_ric.R")

data_gen <- function(n_train,
                     p_feat,
                     par,
                     sigma = 1,
                     eta,
                     seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  qs = par$qs
  qw = par$qw
  cWe = par$cWe
  cWf = par$cWf
  cSe = par$cSe
  cSf = par$cSf
  tau0 = par$tau0
  hetero = par$hetero
  
  ## Covariates
  X_train <- matrix(runif(n_train * p_feat), n_train, p_feat) - 1/2
  
  ## Index sets
  S <- seq_len(qs)
  W <- seq(qs + 1, qs + qw)
  N <- seq(qs + qw + 1, p_feat)
  
  ## Propensity score & treatment
  k <- qs + qw
  scaling_factor <- 1 / sqrt(k)
  
  alpha0 <- 0
  alpha_S <- rep(cSe, qs)
  alpha_W <- rep(cWe, qw)
  
  linpred_e <- alpha0 +
    scaling_factor * (
      X_train[, S, drop = FALSE] %*% alpha_S +
        X_train[, W, drop = FALSE] %*% alpha_W
    )
  
  propensity <- pnorm(linpred_e)
  treatment  <- rbinom(n_train, 1, propensity)
  hist(propensity)
  
  ## Prognostic component
  beta_S <- rep(cSf, qs)
  beta_W <- rep(cWf, qw)
  
  mu <- as.numeric(
    X_train[, S, drop = FALSE] %*% beta_S +
      X_train[, W, drop = FALSE] %*% beta_W
  )
  
  ## Treatment effect
  if (!hetero) {
    tau <- rep(tau0, n_train)
  } else {
    tau_S <- rep(0.5, qs)
    tau <- as.numeric(tau0 + X_train[, S, drop = FALSE] %*% tau_S)
  }
  
  ## Generate log event times
  true_event_times <- mu + (treatment - 0.5) * tau
  uncensored_event_times <- true_event_times +
    rnorm(n_train, 0, sqrt(sigma))
  
  ## Rescale (as in original data_gen)
  sd_un <- sd(uncensored_event_times)
  uncensored_event_times <- uncensored_event_times / sd_un
  true_event_times       <- true_event_times / sd_un
  
  ## Exponential censoring on log scale
  C <- log(rexp(n_train, rate = eta)) + min(uncensored_event_times)
  
  follow_up <- pmin(uncensored_event_times, C)
  status    <- as.numeric(uncensored_event_times <= C)
  
  return(list(
    X_train                = X_train,
    treatment              = treatment,
    propensity             = propensity,
    follow_up              = as.numeric(follow_up),
    status                 = status,
    true_event_times       = as.numeric(true_event_times),
    uncensored_event_times = as.numeric(uncensored_event_times),
    true_cate              = as.vector(tau) / sd_un,
    obs_sigma              = sqrt(sigma) / sd_un,
    sets                   = list(S = S, W = W, N = N),
    params                 = list(
      n_train = n_train,
      p_feat  = p_feat,
      qs      = qs,
      qw      = qw,
      cWe      = cWe,
      cWf      = cWf,
      tau0   = tau0,
      hetero = hetero,
      eta    = eta
    )
  ))
}


## Calibration settings
p_feat      <- 500
n_train     <- 200
N_post      <- 3000
N_burn      <- 2000

sigma       <- 1
eta_hat     <- 0.0137
sim_par <- list(
  qs     = 5,
  qw     = 1,
  cWe     = 5.0,
  cWf     = 1.0,
  cSe    = 1.0,
  cSf    = 5.0,
  tau0   = 1,
  hetero = FALSE
)



# Final simulated dataset
data <- simulate_RIC_data(
  n = n_train,
  p = p_feat,
  par = sim_par,
  eta = eta_hat
)
mean(1-data$status)
mean(data$true_cate)


fit_sBCF <- SurvivalShrinkageBCF(time = data$follow_up,
                                 status = data$status,
                                 X_train = data$X_train,
                                 treatment = data$treatment,
                                 timescale = "log",
                                 propensity = data$propensity,
                                 N_post = N_post,
                                 N_burn = N_burn,
                                 verbose = TRUE
)

rmse_sBCF <- mean((data$true_cate - fit_sBCF$train_predictions_treat)^2)
aBias_sBCF <- mean(abs(data$true_cate - fit_sBCF$train_predictions_treat))

fit_chf <- CausalShrinkageForest(
  y = data$follow_up,
  status = data$status,
  X_train_control = data$X_train,
  X_train_treat = data$X_train,
  treatment_indicator_train = data$treatment,
  outcome_type = "right-censored",
  timescale = "log",
  number_of_trees_control = 200,
  number_of_trees_treat = 200,
  prior_type_control = "horseshoe",
  prior_type_treat = "horseshoe",
  local_hp_control = 0.5 / sqrt(200),
  local_hp_treat = 0.5 / sqrt(200),
  global_hp_control = 0.5 / sqrt(200),
  global_hp_treat = 0.5 / sqrt(200),
  N_post = N_post,
  N_burn = N_burn,
  store_posterior_sample = FALSE,
  verbose = TRUE
)

rmse_chf <- mean((data$true_cate - fit_chf$train_predictions_treat)^2)
aBias_chf <- mean(abs(data$true_cate - fit_chf$train_predictions_treat))

fit_bcf <- SurvivalBCF(time = data$follow_up,
                       status = data$status,
                       X_train = data$X_train,
                       treatment = data$treatment,
                       timescale = "log",
                       propensity = data$propensity,
                       N_post = N_post,
                       N_burn = N_burn,
                       verbose = TRUE
)

rmse_BCF <- mean((data$true_cate - fit_bcf$train_predictions_treat)^2)
aBias_BCF <- mean(abs(data$true_cate - fit_bcf$train_predictions_treat))


results <- data.frame(
  Method = c("Shrinkage BCF", "Causal Horseshoe Forest", "BCF"),
  RMSE   = c(rmse_sBCF, rmse_chf, rmse_BCF),
  AbsBias = c(aBias_sBCF, aBias_chf, aBias_BCF)
)

print(results, row.names = FALSE)


