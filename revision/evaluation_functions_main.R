### Evaluation functions for the simulations ###

library(ShrinkageTrees)

# ============================================================
# Metric helpers
# ============================================================

censored_info <- function(y, status) {
  
  if (requireNamespace("survival", quietly = TRUE)) {
    
    # Fit Gaussian AFT with censoring
    fit <- survival::survreg(survival::Surv(y, status) ~ 1, dist = "gaussian")
    
    mu <- as.numeric(fit$coefficients[1])
    sd <- as.numeric(fit$scale)
    
    # ---- IMPUTE censored values ----
    cens_idx <- which(status == 0)
    imputed_y <- y
    
    if (length(cens_idx) > 0) {
      a <- (y[cens_idx] - mu) / sd
      lambda <- dnorm(a) / (1 - pnorm(a))     # Mills ratio
      imputed_y[cens_idx] <- mu + sd * lambda # Replace censored with E[Y | Y > c]
    }
    
    # min/max after imputation
    min_val <- min(imputed_y, na.rm = TRUE)
    max_val <- max(imputed_y, na.rm = TRUE)
    
  } else {
    
    # Fallback: use only uncensored observations
    min_val <- min(y, na.rm = TRUE)
    max_val <- max(y, na.rm = TRUE)
    
    mu <- mean(y[status == 1], na.rm = TRUE)
    sd <- sd(y[status == 1], na.rm = TRUE)
  }
  
  return(list(
    mu = mu,
    sd = sd,
    min = min_val,
    max = max_val
  ))
}

compute_cate_detection <- function(cate_ci, true_cate, delta) {
  idx <- abs(true_cate) > delta
  if (!any(idx)) return(NA_real_)
  detected <- cate_ci[1, idx] > 0 | cate_ci[2, idx] < 0
  mean(detected)
}

compute_all_metrics <- function(
    cate_samples,
    cate_mean,
    ate_samples,
    train_predictions,
    data,
    alpha = 0.05
) {
  # ---- CATE accuracy ----
  cate_rpehe <- sqrt(mean((cate_mean - data$true_cate)^2))
  cate_abs_bias <- mean(abs(cate_mean - data$true_cate))
  
  # ---- CATE uncertainty ----
  cate_ci <- apply(cate_samples, 2, quantile, probs = c(alpha / 2, 1 - alpha / 2))
  lower <- cate_ci[1, ]
  upper <- cate_ci[2, ]
  
  cate_coverage <- mean(lower <= data$true_cate & upper >= data$true_cate)
  cate_ci_length <- mean(upper - lower)
  
  # ---- Empirical power (CI excludes 0) ----
  empirical_power <- mean(lower > 0 | upper < 0)
  
  # ---- Correct Significant Discovery (CSD) ----
  csd <- mean(
    (lower > 0 | upper < 0) &
      (lower <= data$true_cate & upper >= data$true_cate)
  )
  
  # ---- Winkler score ----
  winkler <- mean(
    (upper - lower) +
      (2 / alpha) * (lower - data$true_cate) * (data$true_cate < lower) +
      (2 / alpha) * (data$true_cate - upper) * (data$true_cate > upper)
  )
  
  # ---- CATE detection power (thresholded) ----
  delta <- 0.25 * sd(data$true_cate)
  cate_detection_power0.25 <- compute_cate_detection(
    cate_ci   = cate_ci,
    true_cate = data$true_cate,
    delta     = delta
  )
  
  delta <- 0.1 * sd(data$true_cate)
  cate_detection_power0.1 <- compute_cate_detection(
    cate_ci   = cate_ci,
    true_cate = data$true_cate,
    delta     = delta
  )
  
  # ---- ATE metrics ----
  ate_estimate <- mean(ate_samples)
  ate_ci <- quantile(ate_samples, probs = c(alpha / 2, 1 - alpha / 2))
  ate_coverage <- ate_ci[1] <= data$true_ate & ate_ci[2] >= data$true_ate
  ate_ci_length <- diff(ate_ci)
  
  # ---- Outcome prediction ----
  rmse <- sqrt(mean((train_predictions - data$true_event_times)^2))
  
  # ---- C-index ----
  if (requireNamespace("survival", quietly = TRUE)) {
    c_index <- survival::concordance(
      survival::Surv(exp(data$follow_up), data$status) ~
        exp(train_predictions)
    )$concordance
  } else {
    c_index <- NA_real_
  }
  
  list(
    ate                   = ate_estimate,
    ate_coverage          = ate_coverage,
    ate_ci_length         = ate_ci_length,
    cate_rpehe            = cate_rpehe,
    cate_abs_bias         = cate_abs_bias,
    cate_coverage         = cate_coverage,
    cate_ci_length        = cate_ci_length,
    empirical_power       = empirical_power,
    csd                   = csd,
    winkler_score         = winkler,
    cate_detection_power0.1  = cate_detection_power0.1,
    cate_detection_power0.25  = cate_detection_power0.25,
    rmse                  = rmse,
    c_index               = c_index
  )
}


evaluate_AFT_BCF_new <- function(
    data,
    use_uncensored = FALSE,
    N_post = N_post,
    N_burn = N_burn,
    ...
) {
  if (!requireNamespace("ShrinkageTrees", quietly = TRUE)) {
    stop("Package 'ShrinkageTrees' is required.")
  }
  
  # --------------------------------------------------
  # Outcome and censoring
  # --------------------------------------------------
  y <- if (use_uncensored) data$uncensored_event_times else data$follow_up
  status <- if (use_uncensored) NULL else data$status
  
  # --------------------------------------------------
  # Fit model
  # --------------------------------------------------
  X_train_control <- cbind(data$propensity, data$X_train)
  
  info <- censored_info(y, status)
  sd_time <- info$sd
  
  fit <- ShrinkageTrees::CausalShrinkageForest(
    y = y,
    status = status,
    outcome_type = "right-censored",
    timescale = "log",
    
    X_train_control = X_train_control,
    X_train_treat   = data$X_train,
    treatment_indicator_train = data$treatment,
    
    prior_type_control = "standard-halfcauchy",
    prior_type_treat = "standard",
    local_hp_treat = 1 / sqrt(50),
    
    power_control = 2,
    base_control  = 0.95,
    power_treat   = 3,
    base_treat    = 0.25,
    
    number_of_trees_control = 200,
    number_of_trees_treat   = 50,
    
    store_posterior_sample = TRUE,
    N_post = N_post,
    N_burn = N_burn,
    verbose = FALSE,
    ...
  )
  
  # --------------------------------------------------
  # Extract posterior quantities
  # --------------------------------------------------
  cate_samples <- fit$train_predictions_sample_treat
  cate_mean    <- fit$train_predictions_treat
  ate_samples  <- rowMeans(cate_samples)
  
  # --------------------------------------------------
  # Compute all metrics (single source of truth)
  # --------------------------------------------------
  metrics <- compute_all_metrics(
    cate_samples      = cate_samples,
    cate_mean         = cate_mean,
    ate_samples       = ate_samples,
    train_predictions = fit$train_predictions,
    data              = data
  )
  
  # --------------------------------------------------
  # Return
  # --------------------------------------------------
  c(
    metrics,
    list(postmean_sigma = mean(fit$sigma))
  )
}


# ============================================================
# AFT-BART
# ============================================================

evaluate_AFT_BART <- function(data, use_uncensored = FALSE,
                              N_post = N_post, N_burn = N_burn, ...) {
  
  y <- if (use_uncensored) data$uncensored_event_times else exp(data$follow_up)
  status <- if (use_uncensored) NULL else data$status
  
  n <- length(y)
  X_train <- cbind(data$treatment, data$X_train)
  X_test <- cbind(rep(c(1, 0), each = n),
                  rbind(data$X_train, data$X_train))
  
  fit <- ShrinkageTrees::SurvivalBART(
    time = y,
    status = status,
    X_train = X_train,
    X_test = X_test,
    timescale = "log",
    N_post = N_post,
    N_burn = N_burn,
    verbose = FALSE,
    ...
  )
  
  cate_samples <- fit$test_predictions_sample[, 1:n] -
    fit$test_predictions_sample[, (n + 1):(2 * n)]
  
  cate_mean  <- colMeans(cate_samples)
  ate_samples <- rowMeans(cate_samples)
  
  metrics <- compute_all_metrics(
    cate_samples      = cate_samples,
    cate_mean         = cate_mean,
    ate_samples       = ate_samples,
    train_predictions = fit$train_predictions,
    data              = data
  )
  
  c(metrics, list(postmean_sigma = mean(fit$sigma)))
}

# ============================================================
# AFT-DART
# ============================================================

evaluate_AFT_DART <- function(data, use_uncensored = FALSE,
                              N_post = N_post, N_burn = N_burn, ...) {
  
  y <- if (use_uncensored) data$uncensored_event_times else exp(data$follow_up)
  status <- if (use_uncensored) NULL else data$status
  
  n <- length(y)
  X_train <- cbind(data$treatment, data$X_train)
  X_test <- cbind(rep(c(1, 0), each = n),
                  rbind(data$X_train, data$X_train))
  
  fit <- ShrinkageTrees::SurvivalDART(
    time = y,
    status = status,
    X_train = X_train,
    X_test = X_test,
    timescale = "log",
    N_post = N_post,
    N_burn = N_burn,
    verbose = FALSE,
    ...
  )
  
  cate_samples <- fit$test_predictions_sample[, 1:n] -
    fit$test_predictions_sample[, (n + 1):(2 * n)]
  
  cate_mean  <- colMeans(cate_samples)
  ate_samples <- rowMeans(cate_samples)
  
  metrics <- compute_all_metrics(
    cate_samples      = cate_samples,
    cate_mean         = cate_mean,
    ate_samples       = ate_samples,
    train_predictions = fit$train_predictions,
    data              = data
  )
  
  c(metrics, list(postmean_sigma = mean(fit$sigma)))
}

# ============================================================
# AFT-BCF
# ============================================================

evaluate_AFT_BCF <- function(data, use_uncensored = FALSE,
                             N_post = N_post, N_burn = N_burn, ...) {
  
  y <- if (use_uncensored) data$uncensored_event_times else data$follow_up
  status <- if (use_uncensored) NULL else data$status
  
  fit <- ShrinkageTrees::SurvivalBCF(
    time = y,
    status = status,
    X_train = data$X_train,
    treatment = data$treatment,
    propensity = data$propensity,
    timescale = "log",
    N_post = N_post,
    N_burn = N_burn,
    verbose = FALSE,
    ...
  )
  
  cate_samples <- fit$train_predictions_sample_treat
  cate_mean   <- fit$train_predictions_treat
  ate_samples <- rowMeans(cate_samples)
  
  metrics <- compute_all_metrics(
    cate_samples      = cate_samples,
    cate_mean         = cate_mean,
    ate_samples       = ate_samples,
    train_predictions = fit$train_predictions,
    data              = data
  )
  
  c(metrics, list(postmean_sigma = mean(fit$sigma)))
}

# ============================================================
# Causal Horseshoe Forest (CHF)
# ============================================================

evaluate_CHF <- function(data, use_uncensored = FALSE,
                         N_post = N_post, N_burn = N_burn,
                         sigma_known = FALSE, ...) {
  
  y <- if (use_uncensored) data$uncensored_event_times else data$follow_up
  status <- if (use_uncensored) NULL else data$status
  sigma <- if (sigma_known) data$obs_sigma else NULL
  
  X_control <- cbind(data$propensity, data$X_train)
  
  fit <- ShrinkageTrees::CausalShrinkageForest(
    y = y,
    status = status,
    outcome_type = "right-censored",
    X_train_control = X_control,
    X_train_treat = data$X_train,
    treatment_indicator_train = data$treatment,
    prior_type_control = "horseshoe",
    prior_type_treat = "horseshoe",
    sigma = sigma,
    N_post = N_post,
    N_burn = N_burn,
    timescale = "log",
    store_posterior_sample = TRUE,
    verbose = FALSE,
    ...
  )
  
  cate_samples <- fit$train_predictions_sample_treat
  cate_mean   <- fit$train_predictions_treat
  ate_samples <- rowMeans(cate_samples)
  
  metrics <- compute_all_metrics(
    cate_samples      = cate_samples,
    cate_mean         = cate_mean,
    ate_samples       = ate_samples,
    train_predictions = fit$train_predictions,
    data              = data
  )
  
  c(metrics, list(postmean_sigma = mean(fit$sigma)))
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
    stop("Package 'survival' is required.")
  }
  
  # --------------------------------------------------
  # Outcome
  # --------------------------------------------------
  y <- if (use_uncensored) data$uncensored_event_times else data$follow_up
  status <- if (use_uncensored) rep(1, length(y)) else data$status
  sigma <- if (sigma_known) data$obs_sigma else NULL
  
  X_control <- cbind(data$propensity, data$X_train)
  X_treat   <- data$X_train
  treatment <- data$treatment
  
  n <- length(y)
  fold_ids <- sample(rep(1:K, length.out = n))
  
  cv_results <- lapply(seq_len(nrow(param_grid)), function(i) {
    k_val <- param_grid$k[i]
    hp <- k_val / sqrt(number_of_trees)
    cindex_fold <- numeric(K)
    
    for (fold in seq_len(K)) {
      test_idx  <- which(fold_ids == fold)
      train_idx <- setdiff(seq_len(n), test_idx)
      
      fit <- ShrinkageTrees::CausalShrinkageForest(
        y = y[train_idx],
        status = status[train_idx],
        outcome_type = "right-censored",
        X_train_control = X_control[train_idx, , drop = FALSE],
        X_train_treat   = X_treat[train_idx, , drop = FALSE],
        X_test_control  = X_control[test_idx, , drop = FALSE],
        X_test_treat    = X_treat[test_idx, , drop = FALSE],
        treatment_indicator_train = treatment[train_idx],
        treatment_indicator_test  = treatment[test_idx],
        number_of_trees_control = number_of_trees,
        number_of_trees_treat   = number_of_trees,
        local_hp_control  = hp,
        global_hp_control = hp,
        local_hp_treat    = hp,
        global_hp_treat   = hp,
        prior_type_control = "horseshoe",
        prior_type_treat   = "horseshoe",
        N_post = N_post,
        N_burn = N_burn,
        sigma = sigma,
        timescale = "log",
        verbose = FALSE
      )
      
      pred <- exp(fit$test_predictions)
      cindex_fold[fold] <- survival::concordance(
        survival::Surv(exp(y[test_idx]), status[test_idx]) ~ pred
      )$concordance
    }
    
    data.frame(k = k_val,
               Mean_CIndex = mean(cindex_fold),
               SD_CIndex   = sd(cindex_fold))
  })
  
  cv_summary <- do.call(rbind, cv_results)
  best_k <- cv_summary$k[which.max(cv_summary$Mean_CIndex)]
  hp <- best_k / sqrt(number_of_trees)
  
  # --------------------------------------------------
  # Refit on full data
  # --------------------------------------------------
  final_fit <- ShrinkageTrees::CausalShrinkageForest(
    y = y,
    status = status,
    outcome_type = "right-censored",
    X_train_control = X_control,
    X_train_treat   = X_treat,
    treatment_indicator_train = treatment,
    number_of_trees_control = number_of_trees,
    number_of_trees_treat   = number_of_trees,
    local_hp_control  = hp,
    global_hp_control = hp,
    local_hp_treat    = hp,
    global_hp_treat   = hp,
    prior_type_control = "horseshoe",
    prior_type_treat   = "horseshoe",
    N_post = N_post,
    N_burn = N_burn,
    sigma = sigma,
    timescale = "log",
    store_posterior_sample = TRUE,
    verbose = FALSE
  )
  
  cate_samples <- final_fit$train_predictions_sample_treat
  cate_mean    <- final_fit$train_predictions_treat
  ate_samples  <- rowMeans(cate_samples)
  
  metrics <- compute_all_metrics(
    cate_samples      = cate_samples,
    cate_mean         = cate_mean,
    ate_samples       = ate_samples,
    train_predictions = final_fit$train_predictions,
    data              = data
  )
  
  c(metrics, list(postmean_sigma = mean(final_fit$sigma), k_cv = best_k))
}


evaluate_AFT_ShrinkageBCF <- function(
    data,
    use_uncensored = FALSE,
    N_post = N_post,
    N_burn = N_burn,
    ...
) {
  if (!requireNamespace("ShrinkageTrees", quietly = TRUE)) {
    stop("Package 'ShrinkageTrees' is required.")
  }
  
  y <- if (use_uncensored) data$uncensored_event_times else data$follow_up
  status <- if (use_uncensored) NULL else data$status
  
  fit <- ShrinkageTrees::SurvivalShrinkageBCF(
    time = y,
    status = status,
    X_train = data$X_train,
    treatment = data$treatment,
    propensity = data$propensity,
    timescale = "log",
    N_post = N_post,
    N_burn = N_burn,
    verbose = FALSE,
    ...
  )
  
  cate_samples <- fit$train_predictions_sample_treat
  cate_mean    <- fit$train_predictions_treat
  ate_samples  <- rowMeans(cate_samples)
  
  metrics <- compute_all_metrics(
    cate_samples      = cate_samples,
    cate_mean         = cate_mean,
    ate_samples       = ate_samples,
    train_predictions = fit$train_predictions,
    data              = data
  )
  
  c(
    metrics,
    list(postmean_sigma = mean(fit$sigma))
  )
}
# ============================================================
# Simulation wrappers (unchanged interface)
# ============================================================

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
  
  res <- list(
    # CHF_CV = evaluate_CHF_CV(
    #   dt,
    #   param_grid = param_grid,
    #   N_post = N_post,
    #   N_burn = N_burn
    # ),
    # CHF0.1 = evaluate_CHF(
    #   dt,
    #   local_hp_control  = 0.1 / sqrt(200),
    #   local_hp_treat    = 0.1 / sqrt(200),
    #   global_hp_control = 0.1 / sqrt(200),
    #   global_hp_treat   = 0.1 / sqrt(200),
    #   N_post = N_post, N_burn = N_burn
    # ),
    CHF_0.2 = evaluate_CHF(
      dt,
      local_hp_control  = 0.2 / sqrt(200),
      local_hp_treat    = 0.2 / sqrt(200),
      global_hp_control = 0.2 / sqrt(200),
      global_hp_treat   = 0.2 / sqrt(200),
      N_post = N_post, N_burn = N_burn
    ),
    CHF0.25 = evaluate_CHF(
      dt,
      local_hp_control  = 0.25 / sqrt(200),
      local_hp_treat    = 0.25 / sqrt(200),
      global_hp_control = 0.25 / sqrt(200),
      global_hp_treat   = 0.25 / sqrt(200),
      N_post = N_post, N_burn = N_burn
    ),
    CHF0.3 = evaluate_CHF(
      dt,
      local_hp_control  = 0.3 / sqrt(200),
      local_hp_treat    = 0.3 / sqrt(200),
      global_hp_control = 0.3 / sqrt(200),
      global_hp_treat   = 0.3 / sqrt(200),
      N_post = N_post, N_burn = N_burn
    ),
    CHF0.5 = evaluate_CHF(
      dt,
      local_hp_control  = 0.5 / sqrt(200),
      local_hp_treat    = 0.5 / sqrt(200),
      global_hp_control = 0.5 / sqrt(200),
      global_hp_treat   = 0.5 / sqrt(200),
      N_post = N_post, N_burn = N_burn
    ),
    AFT_BART = evaluate_AFT_BART(
      dt,
      N_post = N_post,
      N_burn = N_burn
    ),

    AFT_DART = evaluate_AFT_DART(
      dt,
      N_post = N_post,
      N_burn = N_burn
    ),
    AFT_BCF = evaluate_AFT_BCF(
      dt,
      N_post = N_post,
      N_burn = N_burn
    ),
    AFT_ShrinkageBCF = evaluate_AFT_ShrinkageBCF(
      dt,
      N_post = N_post,
      N_burn = N_burn
    )
  )
  
  
  extract <- function(res) {
    data.frame(
      sigma                 = res$postmean_sigma,
      ATE                   = res$ate,
      ATE_bias              = res$ate - dt$true_ate,
      ATE_coverage          = res$ate_coverage,
      ATE_CI_Length         = res$ate_ci_length,
      CATE_RPEHE            = res$cate_rpehe,
      CATE_AbsBias          = res$cate_abs_bias,
      CATE_coverage         = res$cate_coverage,
      CATE_CI_Length        = res$cate_ci_length,
      CATE_Detection_Power0.25  = res$cate_detection_power0.25,
      CATE_Detection_Power0.1  = res$cate_detection_power0.1,
      CATE_winkler  = res$winkler_score,
      CATE_csd  = res$csd,
      CATE_empirical_power  = res$empirical_power,
      RMSE                  = res$rmse,
      C_Index               = res$c_index
    )
  }
  
  out <- do.call(rbind, lapply(res, extract))
  out$Method <- names(res)
  rownames(out) <- NULL
  out
}


# ============================================================
# Parallel wrapper
# ============================================================

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
  results <- foreach::foreach(
    i = 1:M,
    .combine = rbind,
    .packages = c("ShrinkageTrees", "survival"),
    .export = c(
      # generators + metrics
      "data_gen",
      "compute_cate_detection",
      "compute_all_metrics",
      # evaluators
      "evaluate_AFT_BART",
      "evaluate_AFT_DART",
      "evaluate_AFT_BCF",
      "evaluate_AFT_ShrinkageBCF",
      "evaluate_AFT_BCF_new",
      "evaluate_CHF",
      "evaluate_CHF_CV",
      "censored_info",
      # simulation driver
      "run_single_simulation"
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
  
  results
}

# ============================================================
# Over-p wrapper (unchanged interface)
# ============================================================

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
    cat("Running simulations for p =", p_val,
        "cens_scale =", cs_val, "\n")
    
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
  }, p_values, cens_scales)
  
  do.call(rbind, all_results)
}