library(bartcs)
library(foreach)
library(doParallel)
library(grf)

data(ihdp)

# 1. Setup Data
ihdp_cate <- ihdp$mu1 - ihdp$mu0
ihdp_X    <- data.matrix(ihdp[, 6:30]) # Keep as matrix for GRF
ihdp_trt  <- ihdp$treatment
ihdp_y    <- ihdp$y_factual
n         <- length(ihdp_y)

# 2. Compute Propensity Score (Fixing the error)
# We must create a data frame that contains both the predictors and the treatment
# so that glm(treatment ~ .) knows what to do.
df_for_ps <- as.data.frame(ihdp_X)
df_for_ps$treatment <- ihdp_trt

ps_model  <- glm(treatment ~ ., data = df_for_ps, family = binomial(link = "probit"))
ihdp_prop <- predict(ps_model, type = "response")
## ---- helpers: metrics ----
cate_metrics_from_ci <- function(tau_hat, lower, upper, tau_true) {
  rmse <- sqrt(mean((tau_hat - tau_true)^2))
  coverage <- mean(tau_true >= lower & tau_true <= upper)
  int_len <- mean(upper - lower)
  
  data.frame(RMSE = rmse, Coverage = coverage, IntLen = int_len)
}

## ---- GRF fit: returns tau_hat + tau_se ----
fit_cf <- function(p_noise = 0, seed) {
  
  set.seed(seed)
  
  X_train <- cbind(
    ihdp_X,
    matrix(runif(n * p_noise), nrow = n, ncol = p_noise)
  )
  
  cf <- grf::causal_forest(
    X = X_train,
    Y = ihdp_y,
    W = ihdp_trt,
    W.hat = ihdp_prop,
    num.trees = 4000,
    seed = seed
  )
  
  pred <- predict(cf, estimate.variance = TRUE)
  
  list(
    tau_hat = as.numeric(pred$predictions),
    tau_se  = sqrt(as.numeric(pred$variance.estimates))
  )
}

## ---- CHF fit: returns tau_hat + tau_post (for CI) ----
fit_chf <- function(p_noise = 0, k = 0.1, seed,
                    N_post = 2500, N_burn = 2500) {
  
  set.seed(seed)
  
  X_train <- cbind(
    ihdp_X,
    matrix(runif(n * p_noise), nrow = n, ncol = p_noise)
  )
  
  chf <- ShrinkageTrees::CausalShrinkageForest(
    y = ihdp_y,
    X_train_control = cbind(X_train, ihdp_prop),
    X_train_treat   = X_train,
    treatment_indicator_train = ihdp_trt,
    outcome_type = "continuous",
    number_of_trees_control = 200,
    number_of_trees_treat   = 200,
    prior_type_control = "horseshoe",
    prior_type_treat   = "horseshoe",
    local_hp_control  = k/sqrt(200),
    local_hp_treat    = k/sqrt(200),
    global_hp_control = k/sqrt(200),
    global_hp_treat   = k/sqrt(200),
    N_post = N_post,
    N_burn = N_burn,
    store_posterior_sample = TRUE,
    verbose = FALSE
  )
  
  list(
    tau_hat  = as.numeric(chf$train_predictions_treat),
    tau_post = chf$train_predictions_sample_treat
  )
}

## ---- experiment runner ----
run_experiment <- function(p_vals, k_vals, alpha = 0.05, seed = 123) {
  
  z_alpha <- qnorm(1 - alpha / 2)
  out <- list()
  
  for (p in p_vals) {
    cat("Starting scenario: p = ", p, ".\n")
    
    fit <- fit_cf(p_noise = p, seed = seed)
    lower <- fit$tau_hat - z_alpha * fit$tau_se
    upper <- fit$tau_hat + z_alpha * fit$tau_se
    
    met <- cate_metrics_from_ci(
      tau_hat = fit$tau_hat,
      lower = lower,
      upper = upper,
      tau_true = ihdp_cate
    )
    
    out[[length(out) + 1]] <- cbind(
      Method = "Causal Forest",
      p_noise = p,
      k = NA_real_,
      met
    )
    
    for (k in k_vals) {
      
      fit <- fit_chf(p_noise = p, k = k, seed = seed)
      
      ci <- apply(
        fit$tau_post, 2, quantile,
        probs = c(alpha / 2, 1 - alpha / 2),
        na.rm = TRUE
      )
      lower <- ci[1, ]
      upper <- ci[2, ]
      
      met <- cate_metrics_from_ci(
        tau_hat = fit$tau_hat,
        lower = lower,
        upper = upper,
        tau_true = ihdp_cate
      )
      
      out[[length(out) + 1]] <- cbind(
        Method = "Causal Horseshoe Forest",
        p_noise = p,
        k = k,
        met
      )
    }
  }
  
  res <- do.call(rbind, out)
  res[order(res$Method, res$p_noise, res$k), ]
}

## ---- run + print ----
p_vals <- c(0, 500, 1000, 5000)
k_vals <- c(0.1, 0.5, 1, 1.5)

res <- run_experiment(p_vals, k_vals, alpha = 0.05, seed = 123)
print(res, row.names = FALSE)



