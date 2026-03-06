### Evaluation functions for the simulations ###

library(ShrinkageTrees)


compute_cate_detection <- function(cate_ci, true_cate, delta) {
  idx <- abs(true_cate) > delta
  if (!any(idx)) return(NA_real_)
  detected <- cate_ci[1, idx] > 0 | cate_ci[2, idx] < 0
  mean(detected)
}


evaluate_AFT_BART <- function(data, use_uncensored = FALSE, N_post = N_post,
                              N_burn = N_burn, ...) {
  if (!requireNamespace("ShrinkageTrees", quietly = TRUE)) {
    stop("Package 'ShrinkageTrees' is required.")
  }
  
  # Determine outcome and censoring indicator
  if (use_uncensored) {
    y <- data$uncensored_event_times
    status <- NULL
  } else {
    y <- exp(data$follow_up)
    status <- data$status
  }
  
  n <- length(y)
  X_train <- cbind(data$treatment, data$X_train)
  pseudo_treat <- c(rep(1, n), rep(0, n))
  pseudo_cov <- rbind(data$X_train, data$X_train)
  X_test <- cbind(pseudo_treat, pseudo_cov)
  
  # Fit the CHF model
  fit <- ShrinkageTrees::SurvivalBART(
    time = y,
    status = status,
    X_train = X_train,
    X_test = X_test,
    timescale = "log",
    verbose = FALSE,
    N_post = N_post,
    N_burn = N_burn,
    ...
  )
  
  # Posterior samples of τ(x)
  cate_samples <- fit$test_predictions_sample[,c(1:n)] - fit$test_predictions_sample[,c((n+1):(2*n))]
  
  # Pointwise quantities
  cate_mean <- colMeans(cate_samples)
  cate_rpehe <- sqrt(mean((cate_mean - data$true_cate)^2))
  
  # Per-observation posterior intervals for τ(x)
  cate_ci <- apply(cate_samples, 2, quantile, probs = c(0.025, 0.975))
  cate_coverage <- mean(cate_ci[1,] <= data$true_cate & cate_ci[2,] >= data$true_cate)
  cate_ci_length <- mean(cate_ci[2,] - cate_ci[1,])
  
  # ATE estimate and uncertainty
  ate_samples <- rowMeans(cate_samples)
  ate_estimate <- mean(ate_samples)
  ate_ci <- quantile(ate_samples, probs = c(0.025, 0.975))
  ate_coverage <- unname(ate_ci[1] <= data$true_ate & ate_ci[2] >= data$true_ate)
  ate_ci_length <- unname(diff(ate_ci))
  
  delta <- 0.25 * sd(data$true_cate)
  cate_detection_power <- compute_cate_detection(
    cate_ci   = cate_ci,
    true_cate = data$true_cate,
    delta     = delta
  )
  
  # RMSE for total outcome prediction
  rmse <- sqrt(mean((fit$train_predictions - data$true_event_times)^2))
  
  # Compute the C-index
  if (requireNamespace("survival", quietly = TRUE)) {
    c_index <- survival::concordance(
      survival::Surv(exp(data$follow_up), data$status) ~ exp(fit$train_predictions)
    )$concordance
  } else {
    c_index <- NA
  }
  
  return(list(
    postmean_sigma = mean(fit$sigma),
    ate = ate_estimate,
    ate_ci = ate_ci,
    ate_coverage = ate_coverage,
    ate_ci_length = ate_ci_length,
    cate_rpehe = cate_rpehe,
    cate_coverage = cate_coverage,
    cate_ci_length = cate_ci_length,
    rmse = rmse,
    cate_detection_power  = cate_detection_power,
    c_index = c_index
  ))
}

evaluate_AFT_DART <- function(data, use_uncensored = FALSE, N_post = N_post,
                              N_burn = N_burn, ...) {
  if (!requireNamespace("ShrinkageTrees", quietly = TRUE)) {
    stop("Package 'ShrinkageTrees' is required.")
  }
  
  # Determine outcome and censoring indicator
  if (use_uncensored) {
    y <- data$uncensored_event_times
    status <- NULL
  } else {
    y <- exp(data$follow_up)
    status <- data$status
  }
  
  n <- length(y)
  X_train <- cbind(data$treatment, data$X_train)
  pseudo_treat <- c(rep(1, n), rep(0, n))
  pseudo_cov <- rbind(data$X_train, data$X_train)
  X_test <- cbind(pseudo_treat, pseudo_cov)
  
  # Fit the CHF model
  fit <- ShrinkageTrees::SurvivalDART(
    time = y,
    status = status,
    X_train = X_train,
    X_test = X_test,
    timescale = "log",
    verbose = FALSE,
    N_post = N_post,
    N_burn = N_burn,
    ...
  )
  
  # Posterior samples of τ(x)
  cate_samples <- fit$test_predictions_sample[,c(1:n)] - fit$test_predictions_sample[,c((n+1):(2*n))]
  
  # Pointwise quantities
  cate_mean <- colMeans(cate_samples)
  cate_rpehe <- sqrt(mean((cate_mean - data$true_cate)^2))
  
  # Per-observation posterior intervals for τ(x)
  cate_ci <- apply(cate_samples, 2, quantile, probs = c(0.025, 0.975))
  cate_coverage <- mean(cate_ci[1,] <= data$true_cate & cate_ci[2,] >= data$true_cate)
  cate_ci_length <- mean(cate_ci[2,] - cate_ci[1,])
  
  # ATE estimate and uncertainty
  ate_samples <- rowMeans(cate_samples)
  ate_estimate <- mean(ate_samples)
  ate_ci <- quantile(ate_samples, probs = c(0.025, 0.975))
  ate_coverage <- unname(ate_ci[1] <= data$true_ate & ate_ci[2] >= data$true_ate)
  ate_ci_length <- unname(diff(ate_ci))
  
  delta <- 0.25 * sd(data$true_cate)
  cate_detection_power <- compute_cate_detection(
    cate_ci   = cate_ci,
    true_cate = data$true_cate,
    delta     = delta
  )
  
  # RMSE for total outcome prediction
  rmse <- sqrt(mean((fit$train_predictions - data$true_event_times)^2))
  
  # Compute the C-index
  if (requireNamespace("survival", quietly = TRUE)) {
    c_index <- survival::concordance(
      survival::Surv(exp(data$follow_up), data$status) ~ exp(fit$train_predictions)
    )$concordance
  } else {
    c_index <- NA
  }
  
  return(list(
    postmean_sigma = mean(fit$sigma),
    ate = ate_estimate,
    ate_ci = ate_ci,
    ate_coverage = ate_coverage,
    ate_ci_length = ate_ci_length,
    cate_rpehe = cate_rpehe,
    cate_coverage = cate_coverage,
    cate_ci_length = cate_ci_length,
    rmse = rmse,
    cate_detection_power  = cate_detection_power,
    c_index = c_index
  ))
}


evaluate_AFT_BCF <- function(data, use_uncensored = FALSE, N_post = N_post,
                             N_burn = N_burn, ...) {
  if (!requireNamespace("ShrinkageTrees", quietly = TRUE)) {
    stop("Package 'ShrinkageTrees' is required.")
  }
  
  # Determine outcome and censoring indicator
  if (use_uncensored) {
    y <- data$uncensored_event_times
    status <- NULL
  } else {
    y <- data$follow_up
    status <- data$status
  }
  
  # Fit the CHF model
  fit <- ShrinkageTrees::SurvivalBCF(
    time = y,
    status = status,
    X_train = data$X_train,
    treatment = data$treatment,
    propensity = data$propensity,
    timescale = "log",
    verbose = FALSE,
    N_post = N_post,
    N_burn = N_burn,
    ...
  )
  
  # Posterior samples of τ(x)
  cate_samples <- fit$train_predictions_sample_treat
  
  # Pointwise quantities
  cate_mean <- fit$train_predictions_treat
  cate_rpehe <- sqrt(mean((cate_mean - data$true_cate)^2))
  
  # Per-observation posterior intervals for τ(x)
  cate_ci <- apply(cate_samples, 2, quantile, probs = c(0.025, 0.975))
  cate_coverage <- mean(cate_ci[1,] <= data$true_cate & cate_ci[2,] >= data$true_cate)
  cate_ci_length <- mean(cate_ci[2,] - cate_ci[1,])
  
  # ATE estimate and uncertainty
  ate_samples <- rowMeans(cate_samples)
  ate_estimate <- mean(ate_samples)
  ate_ci <- quantile(ate_samples, probs = c(0.025, 0.975))
  ate_coverage <- unname(ate_ci[1] <= data$true_ate & ate_ci[2] >= data$true_ate)
  ate_ci_length <- unname(diff(ate_ci))
  
  delta <- 0.25 * sd(data$true_cate)
  cate_detection_power <- compute_cate_detection(
    cate_ci   = cate_ci,
    true_cate = data$true_cate,
    delta     = delta
  )
  
  # RMSE for total outcome prediction
  rmse <- sqrt(mean((fit$train_predictions - data$true_event_times)^2))
  
  # Compute the C-index
  if (requireNamespace("survival", quietly = TRUE)) {
    c_index <- survival::concordance(
      survival::Surv(exp(data$follow_up), data$status) ~ exp(fit$train_predictions)
    )$concordance
  } else {
    c_index <- NA
  }
  
  return(list(
    postmean_sigma = mean(fit$sigma),
    ate = ate_estimate,
    ate_ci = ate_ci,
    ate_coverage = ate_coverage,
    ate_ci_length = ate_ci_length,
    cate_rpehe = cate_rpehe,
    cate_coverage = cate_coverage,
    cate_ci_length = cate_ci_length,
    rmse = rmse,
    cate_detection_power  = cate_detection_power,
    c_index = c_index
  ))
}

evaluate_AFT_ShrinkageBCF <- function(data, use_uncensored = FALSE, 
                                      N_post = N_post, N_burn = N_burn, ...) {
  if (!requireNamespace("ShrinkageTrees", quietly = TRUE)) {
    stop("Package 'ShrinkageTrees' is required.")
  }
  
  # Determine outcome and censoring indicator
  if (use_uncensored) {
    y <- data$uncensored_event_times
    status <- NULL
  } else {
    y <- data$follow_up
    status <- data$status
  }
  
  # Fit the CHF model
  fit <- ShrinkageTrees::SurvivalShrinkageBCF(
    time = y,
    status = status,
    X_train = data$X_train,
    treatment = data$treatment,
    propensity = data$propensity,
    timescale = "log",
    verbose = FALSE,
    N_post = N_post,
    N_burn = N_burn,
    ...
  )
  
  # Posterior samples of τ(x)
  cate_samples <- fit$train_predictions_sample_treat
  
  # Pointwise quantities
  cate_mean <- fit$train_predictions_treat
  cate_rpehe <- sqrt(mean((cate_mean - data$true_cate)^2))
  
  # Per-observation posterior intervals for τ(x)
  cate_ci <- apply(cate_samples, 2, quantile, probs = c(0.025, 0.975))
  cate_coverage <- mean(cate_ci[1,] <= data$true_cate & cate_ci[2,] >= data$true_cate)
  cate_ci_length <- mean(cate_ci[2,] - cate_ci[1,])
  
  # ATE estimate and uncertainty
  ate_samples <- rowMeans(cate_samples)
  ate_estimate <- mean(ate_samples)
  ate_ci <- quantile(ate_samples, probs = c(0.025, 0.975))
  ate_coverage <- unname(ate_ci[1] <= data$true_ate & ate_ci[2] >= data$true_ate)
  ate_ci_length <- unname(diff(ate_ci))
  
  delta <- 0.25 * sd(data$true_cate)
  cate_detection_power <- compute_cate_detection(
    cate_ci   = cate_ci,
    true_cate = data$true_cate,
    delta     = delta
  )
  
  # RMSE for total outcome prediction
  rmse <- sqrt(mean((fit$train_predictions - data$true_event_times)^2))
  
  # Compute the C-index
  if (requireNamespace("survival", quietly = TRUE)) {
    c_index <- survival::concordance(
      survival::Surv(exp(data$follow_up), data$status) ~ exp(fit$train_predictions)
    )$concordance
  } else {
    c_index <- NA
  }
  
  return(list(
    postmean_sigma = mean(fit$sigma),
    ate = ate_estimate,
    ate_ci = ate_ci,
    ate_coverage = ate_coverage,
    ate_ci_length = ate_ci_length,
    cate_rpehe = cate_rpehe,
    cate_coverage = cate_coverage,
    cate_ci_length = cate_ci_length,
    rmse = rmse,
    cate_detection_power  = cate_detection_power,
    c_index = c_index
  ))
}


evaluate_CHF <- function(data, use_uncensored = FALSE, N_post = N_post,
                         N_burn = N_burn, sigma_known = FALSE, ...) {
  if (!requireNamespace("ShrinkageTrees", quietly = TRUE)) {
    stop("Package 'ShrinkageTrees' is required.")
  }
  
  # Determine outcome and censoring indicator
  if (use_uncensored) {
    y <- data$uncensored_event_times
    status <- NULL
  } else {
    y <- data$follow_up
    status <- data$status
  }
  
  if (sigma_known) {
    sigma <- data$obs_sigma
  } else {
    sigma <- NULL
  }
  
  X1 <- cbind(data$propensity, data$X_train)
  X2 <- matrix(colMeans(X1), nrow = 1)
  
  # Fit the CHF model
  fit <- ShrinkageTrees::CausalShrinkageForest(
    y = y,
    status = status,
    outcome_type = "right-censored",
    X_train_control = X1,
    X_train_treat = data$X_train,
    treatment_indicator_train = data$treatment,
    store_posterior_sample = TRUE,
    timescale = "log",
    verbose = FALSE,
    prior_type_control = "horseshoe",
    prior_type_treat = "horseshoe",
    sigma = sigma,
    N_post = N_post,
    N_burn = N_burn,
    ...
  )
  
  # Posterior samples of τ(x)
  cate_samples <- fit$train_predictions_sample_treat
  
  # Pointwise quantities
  cate_mean <- fit$train_predictions_treat
  cate_rpehe <- sqrt(mean((cate_mean - data$true_cate)^2))
  
  # Per-observation posterior intervals for τ(x)
  cate_ci <- apply(cate_samples, 2, quantile, probs = c(0.025, 0.975))
  cate_coverage <- mean(cate_ci[1,] <= data$true_cate & cate_ci[2,] >= data$true_cate)
  cate_ci_length <- mean(cate_ci[2,] - cate_ci[1,])
  
  # ATE estimate and uncertainty
  ate_samples <- rowMeans(cate_samples)
  ate_estimate <- mean(ate_samples)
  ate_ci <- quantile(ate_samples, probs = c(0.025, 0.975))
  ate_coverage <- unname(ate_ci[1] <= data$true_ate & ate_ci[2] >= data$true_ate)
  ate_ci_length <- unname(diff(ate_ci))
  
  delta <- 0.25 * sd(data$true_cate)
  cate_detection_power <- compute_cate_detection(
    cate_ci   = cate_ci,
    true_cate = data$true_cate,
    delta     = delta
  )
  
  # RMSE for total outcome prediction
  rmse <- sqrt(mean((fit$train_predictions - data$true_event_times)^2))
  
  # Compute the C-index
  if (requireNamespace("survival", quietly = TRUE)) {
    c_index <- survival::concordance(
      survival::Surv(exp(data$follow_up), data$status) ~ exp(fit$train_predictions)
    )$concordance
  } else {
    c_index <- NA
  }
  
  return(list(
    postmean_sigma = mean(fit$sigma),
    ate = ate_estimate,
    ate_ci = ate_ci,
    ate_coverage = ate_coverage,
    ate_ci_length = ate_ci_length,
    cate_rpehe = cate_rpehe,
    cate_coverage = cate_coverage,
    cate_ci_length = cate_ci_length,
    rmse = rmse,
    cate_detection_power  = cate_detection_power,
    c_index = c_index
  ))
}

evaluate_CHF_CV <- function(data,
                            param_grid,
                            K = 5,
                            number_of_trees = 200,
                            N_post = 1000,
                            N_burn = 1000,
                            use_uncensored = FALSE,
                            sigma_known = FALSE,
                            seed = 1) {
  set.seed(seed)
  
  if (!requireNamespace("ShrinkageTrees", quietly = TRUE)) {
    stop("Package 'ShrinkageTrees' is required.")
  }
  
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package 'survival' is required for computing C-index.")
  }
  
  # Determine outcome
  if (use_uncensored) {
    y_full <- data$uncensored_event_times
    status_full <- rep(1, length(data$status))
  } else {
    y_full <- data$follow_up
    status_full <- data$status
  }
  
  if (sigma_known) {
    sigma <- data$obs_sigma
  } else {
    sigma <- NULL
  }
  
  X_control <- cbind(data$propensity, data$X_train)
  X_treat <- cbind(data$propensity, data$X_train)
  treatment <- data$treatment
  true_event_times <- data$true_event_times
  
  n <- length(y_full)
  fold_ids <- sample(rep(1:K, length.out = n))
  
  results_list <- list()
  
  for (i in seq_len(nrow(param_grid))) {
    k_val <- param_grid$k[i]
    local_hp  <- k_val / sqrt(number_of_trees)
    global_hp <- k_val / sqrt(number_of_trees)
    
    cindex_fold <- numeric(K)
    
    for (fold in 1:K) {
      test_idx <- which(fold_ids == fold)
      train_idx <- setdiff(seq_len(n), test_idx)
      
      y_train <- y_full[train_idx]
      status_train <- status_full[train_idx]
      X_train_control <- X_control[train_idx, , drop = FALSE]
      X_train_treat <- X_treat[train_idx, , drop = FALSE]
      treat_train <- treatment[train_idx]
      
      surv_time_test <- exp(y_full[test_idx])
      status_test <- status_full[test_idx]
      X_test_control <- X_control[test_idx, , drop = FALSE]
      X_test_treat <- X_treat[test_idx, , drop = FALSE]
      treat_test <- treatment[test_idx]
      
      fit <- ShrinkageTrees::CausalShrinkageForest(
        y = y_train,
        status = status_train,
        X_train_control = X_train_control,
        X_train_treat = X_train_treat,
        X_test_control = X_test_control,
        X_test_treat = X_test_treat,
        outcome_type = "right-censored",
        treatment_indicator_train = treat_train,
        treatment_indicator_test = treat_test,
        number_of_trees_control = number_of_trees,
        number_of_trees_treat = number_of_trees,
        local_hp_control = local_hp,
        global_hp_control = global_hp,
        local_hp_treat = local_hp,
        global_hp_treat = global_hp,
        prior_type_control = "horseshoe",
        prior_type_treat = "horseshoe",
        N_post = N_post,
        N_burn = N_burn,
        sigma = sigma,
        verbose = FALSE,
        timescale = "log",
        store_posterior_sample = FALSE
      )
      
      pred <- exp(fit$test_predictions)
      
      c_index <- survival::concordance(
        survival::Surv(surv_time_test, status_test) ~ pred
      )$concordance
      
      cindex_fold[fold] <- c_index
    }
    
    results_list[[i]] <- data.frame(
      k = k_val,
      Mean_CIndex = mean(cindex_fold),
      SD_CIndex = sd(cindex_fold)
    )
  }
  
  cv_summary <- do.call(rbind, results_list)
  
  # Select params with highest mean C-index
  best_row <- cv_summary[which.max(cv_summary$Mean_CIndex), ]
  
  local_hp  <- best_row$k / sqrt(number_of_trees)
  global_hp <- best_row$k / sqrt(number_of_trees)
  
  # Refit on full data using best params
  final_fit <- ShrinkageTrees::CausalShrinkageForest(
    y = y_full,
    status = status_full,
    outcome_type = "right-censored",
    X_train_control = X_control,
    X_train_treat = X_treat,
    treatment_indicator_train = treatment,
    number_of_trees_control = number_of_trees,
    number_of_trees_treat = number_of_trees,
    local_hp_control = local_hp,
    global_hp_control = global_hp,
    local_hp_treat = local_hp,
    global_hp_treat = global_hp,
    prior_type_control = "horseshoe",
    prior_type_treat = "horseshoe",
    N_post = N_post,
    N_burn = N_burn,
    sigma = sigma,
    verbose = FALSE,
    timescale = "log",
    store_posterior_sample = TRUE,
  )
  
  # Posterior samples of τ(x)
  cate_samples <- final_fit$train_predictions_sample_treat
  
  # Pointwise quantities
  cate_mean <- final_fit$train_predictions_treat
  cate_rpehe <- sqrt(mean((cate_mean - data$true_cate)^2))
  
  # Per-observation posterior intervals for τ(x)
  cate_ci <- apply(cate_samples, 2, quantile, probs = c(0.025, 0.975))
  cate_coverage <- mean(cate_ci[1,] <= data$true_cate & cate_ci[2,] >= data$true_cate)
  cate_ci_length <- mean(cate_ci[2,] - cate_ci[1,])
  
  # ATE estimate and uncertainty
  ate_samples <- rowMeans(cate_samples)
  ate_estimate <- mean(ate_samples)
  ate_ci <- quantile(ate_samples, probs = c(0.025, 0.975))
  ate_coverage <- unname(ate_ci[1] <= data$true_ate & ate_ci[2] >= data$true_ate)
  ate_ci_length <- unname(diff(ate_ci))
  
  delta <- 0.25 * sd(data$true_cate)
  cate_detection_power <- compute_cate_detection(
    cate_ci   = cate_ci,
    true_cate = data$true_cate,
    delta     = delta
  )
  
  # RMSE for total outcome prediction
  rmse <- sqrt(mean((final_fit$train_predictions - data$true_event_times)^2))
  
  # Compute the C-index
  if (requireNamespace("survival", quietly = TRUE)) {
    c_index <- survival::concordance(
      survival::Surv(exp(data$follow_up), data$status) ~ exp(final_fit$train_predictions)
    )$concordance
  } else {
    c_index <- NA
  }
  
  return(list(
    postmean_sigma = mean(final_fit$sigma),
    ate = ate_estimate,
    ate_ci = ate_ci,
    ate_coverage = ate_coverage,
    ate_ci_length = ate_ci_length,
    cate_rpehe = cate_rpehe,
    cate_coverage = cate_coverage,
    cate_ci_length = cate_ci_length,
    rmse = rmse,
    cate_detection_power  = cate_detection_power,
    c_index = c_index
  ))
}


run_single_simulation <- function(n_train,
                                  p_feat,
                                  sigma,
                                  s_prog,
                                  s_treat,
                                  cens_scale,
                                  k,
                                  linear,
                                  param_grid,
                                  N_post = 2000,
                                  N_burn = 2000,
                                  seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  dt <- data_gen(n_train, p_feat, sigma, s_prog, s_treat, cens_scale, linear)
  
  res_CHF <- evaluate_CHF(
    dt,
    number_of_trees_control = 200,
    number_of_trees_treat = 200,
    local_hp_control = k / sqrt(200),
    local_hp_treat = k / sqrt(200),
    global_hp_control = k / sqrt(200),
    global_hp_treat = k / sqrt(200),
    N_post = N_post,
    N_burn = N_burn
  )
  
  res_CHF_CV <- evaluate_CHF_CV(
    dt,
    param_grid = param_grid,
    K = 5,
    number_of_trees = 200,
    N_post = N_post,
    N_burn = N_burn
  )
  
  res_CHF_sigma_known <- evaluate_CHF(
    dt,
    number_of_trees_control = 200,
    number_of_trees_treat = 200,
    local_hp_control = k / sqrt(200),
    local_hp_treat = k / sqrt(200),
    global_hp_control = k / sqrt(200),
    global_hp_treat = k / sqrt(200),
    N_post = N_post,
    N_burn = N_burn,
    sigma_known = TRUE
  )
  
  res_CHF_CV_sigma_known <- evaluate_CHF_CV(
    dt,
    param_grid = param_grid,
    K = 5,
    number_of_trees = 200,
    N_post = N_post,
    N_burn = N_burn,
    sigma_known = TRUE
  )
  
  res_AFT_BART <- evaluate_AFT_BART(
    dt,
    N_post = N_post,
    N_burn = N_burn
  )
  
  res_AFT_DART <- evaluate_AFT_DART(
    dt,
    N_post = N_post,
    N_burn = N_burn
  )
  
  res_AFT_BCF <- evaluate_AFT_BCF(
    dt,
    N_post = N_post,
    N_burn = N_burn
  )
  
  res_AFT_S_BCF <- evaluate_AFT_ShrinkageBCF(
    dt,
    N_post = N_post,
    N_burn = N_burn
  )
  
  all_methods <- list(
    CHF            = res_CHF,
    CHF_CV         = res_CHF_CV,
    res_CHF_sigma_known = res_CHF_sigma_known,
    res_CHF_CV_sigma_known = res_CHF_CV_sigma_known,
    AFT_BART       = res_AFT_BART,
    AFT_DART       = res_AFT_DART,
    AFT_BCF        = res_AFT_BCF,
    AFT_S_BCF      = res_AFT_S_BCF
  )
  
  extract <- function(res) {
    data.frame(
      sigma          = res$postmean_sigma,
      ATE            = res$ate,
      ATE_bias       = res$ate - dt$true_ate,
      ATE_coverage   = res$ate_coverage,
      ATE_CI_Length  = res$ate_ci_length,
      CATE_RPEHE     = res$cate_rpehe,
      CATE_coverage  = res$cate_coverage,
      CATE_CI_Length = res$cate_ci_length,
      CATE_Detection_Power  = res$cate_detection_power,
      RMSE           = res$rmse,
      C_Index        = res$c_index
    )
  }
  
  result_df <- do.call(rbind, lapply(all_methods, extract))
  result_df$Method <- names(all_methods)
  rownames(result_df) <- NULL
  result_df <- result_df[, c("Method", setdiff(names(result_df), "Method"))]
  
  return(result_df)
}

run_parallel_simulations <- function(M,
                                     n_train,
                                     p_feat,
                                     sigma,
                                     s_prog,
                                     s_treat,
                                     cens_scale,
                                     k,
                                     linear,
                                     param_grid,
                                     N_post = 2000,
                                     N_burn = 2000) {
  
  results <- foreach(
    i = 1:M,
    .combine = rbind,
    .packages = c("ShrinkageTrees", "survival"),
    .export = c(
      # core components
      "run_single_simulation",
      "data_gen",
      
      # Individual methods
      "evaluate_CHF",
      "evaluate_CHF_CV",
      "evaluate_AFT_BART",
      "evaluate_AFT_DART",
      "evaluate_AFT_BCF",
      "evaluate_AFT_ShrinkageBCF"
    )
  ) %dopar% {
    
    run_single_simulation(
      n_train = n_train,
      p_feat = p_feat,
      sigma = sigma,
      s_prog = s_prog,
      s_treat = s_treat,
      cens_scale = cens_scale,
      k = k,
      linear = linear,
      param_grid = param_grid,
      N_post = N_post,
      N_burn = N_burn,
      seed = i
    )
  }
  
  return(results)
}

run_simulations_over_p <- function(
    p_values,
    cens_scales,
    M,
    n_train,
    sigma,
    s_prog,
    s_treat,
    k,
    linear,
    param_grid,
    N_post = 2000,
    N_burn = 2000 
) {
  
  if (length(p_values) != length(cens_scales)) {
    stop("p_values and cens_scales must be the same length.")
  }
  
  all_results <- Map(function(p_val, cs_val) {
    
    cat("Running simulations for p =", p_val, "cens_scale =", cs_val, "\n")
    
    res <- run_parallel_simulations(
      M = M,
      n_train = n_train,
      p_feat = p_val,
      sigma = sigma,
      s_prog = s_prog,
      s_treat = s_treat,
      cens_scale = cs_val,
      k = k,
      linear = linear,
      param_grid = param_grid,
      N_post = N_post,
      N_burn = N_burn
    )
    
    res$p_feat <- p_val
    res$cens_scale <- cs_val
    res
  },
  p_values, cens_scales)
  
  do.call(rbind, all_results)
}

