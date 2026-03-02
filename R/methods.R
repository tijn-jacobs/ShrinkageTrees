#' @export
predict.ShrinkageTrees <- function(object, newdata, level = 0.95, ...) {

  if (!is.matrix(newdata)) newdata <- as.matrix(newdata)
  if (ncol(newdata) != object$data_info$p_features)
    stop("newdata has ", ncol(newdata), " columns but the model expects ",
         object$data_info$p_features, ".")

  n_new   <- nrow(newdata)
  X_new   <- as.numeric(t(newdata))
  n_train <- object$data_info$n_train
  p       <- object$data_info$p_features
  pre     <- object$preprocess
  args    <- object$args
  alpha   <- (1 - level) / 2

  if (object$outcome_type == "binary") {

    fit <- probitHorseTrees_cpp(
      nSEXP                    = n_train,
      pSEXP                    = p,
      n_testSEXP               = n_new,
      X_trainSEXP              = object$data$X_train,
      ySEXP                    = object$data$y_train,
      X_testSEXP               = X_new,
      number_of_treesSEXP      = object$mcmc$number_of_trees,
      N_postSEXP               = object$mcmc$N_post,
      N_burnSEXP               = object$mcmc$N_burn,
      delayed_proposalSEXP     = args$delayed_proposal,
      powerSEXP                = args$power,
      baseSEXP                 = args$base,
      p_growSEXP               = args$p_grow,
      p_pruneSEXP              = args$p_prune,
      omegaSEXP                = args$omega,
      latent_thresholdSEXP     = pre$latent_threshold,
      param1SEXP               = args$param1,
      param2SEXP               = args$param2,
      prior_typeSEXP           = object$prior$prior_type_cpp,
      reversibleSEXP           = args$reversible,
      store_posterior_sampleSEXP = TRUE,
      verboseSEXP              = FALSE
    )

    out <- list(
      mean         = pnorm(fit$test_predictions),
      lower        = pnorm(apply(fit$test_predictions_sample, 2, quantile, alpha)),
      upper        = pnorm(apply(fit$test_predictions_sample, 2, quantile, 1 - alpha)),
      n            = n_new,
      level        = level,
      outcome_type = object$outcome_type,
      timescale    = object$timescale
    )
    class(out) <- "ShrinkageTreesPrediction"
    return(out)

  } else {

    y_std    <- (object$data$y_train - pre$y_mean) / pre$sigma_hat
    survival <- object$outcome_type == "right-censored"
    status   <- if (survival) object$data$status_train else rep(1L, n_train)

    fit <- HorseTrees_cpp(
      nSEXP                    = n_train,
      pSEXP                    = p,
      n_testSEXP               = n_new,
      X_trainSEXP              = object$data$X_train,
      ySEXP                    = y_std,
      status_indicatorSEXP     = status,
      is_survivalSEXP          = survival,
      X_testSEXP               = X_new,
      number_of_treesSEXP      = object$mcmc$number_of_trees,
      N_postSEXP               = object$mcmc$N_post,
      N_burnSEXP               = object$mcmc$N_burn,
      delayed_proposalSEXP     = args$delayed_proposal,
      powerSEXP                = args$power,
      baseSEXP                 = args$base,
      p_growSEXP               = args$p_grow,
      p_pruneSEXP              = args$p_prune,
      nuSEXP                   = args$nu,
      lambdaSEXP               = args$lambda,
      dirichlet_boolSEXP       = isTRUE(object$prior$dirichlet),
      a_dirichletSEXP          = args$a_dirichlet,
      b_dirichletSEXP          = args$b_dirichlet,
      rho_dirichletSEXP        = args$rho_dirichlet,
      sigmaSEXP                = pre$sigma_hat,
      sigma_knownSEXP          = TRUE,
      omegaSEXP                = args$omega,
      param1SEXP               = args$param1,
      param2SEXP               = args$param2,
      prior_typeSEXP           = object$prior$prior_type_cpp,
      reversibleSEXP           = args$reversible,
      store_parametersSEXP     = FALSE,
      store_posterior_sampleSEXP = TRUE,
      verboseSEXP              = FALSE
    )

    # De-standardize
    if (survival && object$timescale == "time") {
      mean_pred   <- exp(fit$test_predictions       * pre$sigma_hat + pre$y_mean)
      sample_pred <- exp(fit$test_predictions_sample * pre$sigma_hat + pre$y_mean)
    } else {
      mean_pred   <- fit$test_predictions       * pre$sigma_hat + pre$y_mean
      sample_pred <- fit$test_predictions_sample * pre$sigma_hat + pre$y_mean
    }

    out <- list(
      mean         = mean_pred,
      lower        = apply(sample_pred, 2, quantile, alpha),
      upper        = apply(sample_pred, 2, quantile, 1 - alpha),
      n            = n_new,
      level        = level,
      outcome_type = object$outcome_type,
      timescale    = object$timescale
    )
    class(out) <- "ShrinkageTreesPrediction"
    out
  }
}

#' @export
print.ShrinkageTreesPrediction <- function(x, n_head = 6, digits = 3, ...) {

  ci_pct <- paste0(round(x$level * 100), "%")
  scale_label <- switch(
    x$outcome_type,
    "binary"        = "probability",
    "right-censored" = if (x$timescale == "time") "survival time" else "log survival time",
    "fitted value"
  )

  cat("\n")
  cat("ShrinkageTrees predictions\n")
  cat("--------------------------\n")
  lbl <- function(s) sprintf("%-22s", s)
  cat(lbl("Observations:"),      x$n,         "\n", sep = "")
  cat(lbl("Credible interval:"), ci_pct,      "\n", sep = "")
  cat(lbl("Scale:"),             scale_label, "\n", sep = "")

  n_show <- min(n_head, x$n)
  cat("\n")
  cat(sprintf("  %5s  %8s  %8s  %8s\n", "", "mean", "lower", "upper"))
  cat(sprintf("  %5s  %8s  %8s  %8s\n", "", strrep("-", 8), strrep("-", 8), strrep("-", 8)))
  for (i in seq_len(n_show)) {
    cat(sprintf("  [%3d]  %8.*f  %8.*f  %8.*f\n",
                i, digits, x$mean[i], digits, x$lower[i], digits, x$upper[i]))
  }
  if (x$n > n_head)
    cat("  ... (", x$n - n_head, " more)\n", sep = "")

  cat("\n")
  invisible(x)
}

#' @export
summary.ShrinkageTreesPrediction <- function(object, ...) {
  out <- list(
    n            = object$n,
    level        = object$level,
    outcome_type = object$outcome_type,
    timescale    = object$timescale,
    mean         = quantile(object$mean,  c(0, 0.25, 0.5, 1)),
    lower        = quantile(object$lower, c(0, 0.25, 0.5, 1)),
    upper        = quantile(object$upper, c(0, 0.25, 0.5, 1))
  )
  class(out) <- "summary.ShrinkageTreesPrediction"
  out
}

#' @export
print.summary.ShrinkageTreesPrediction <- function(x, digits = 3, ...) {

  ci_pct <- paste0(round(x$level * 100), "%")

  cat("\n")
  cat("ShrinkageTrees prediction summary\n")
  cat("----------------------------------\n")
  lbl <- function(s) sprintf("%-22s", s)
  cat(lbl("Observations:"),      x$n,    "\n", sep = "")
  cat(lbl("Credible interval:"), ci_pct, "\n", sep = "")

  cat("\n")
  cat(sprintf("  %-8s  %8s  %8s  %8s\n", "", "mean", "lower", "upper"))
  cat(sprintf("  %-8s  %8s  %8s  %8s\n", "", strrep("-", 8), strrep("-", 8), strrep("-", 8)))
  nms <- c("Min.", "Q1", "Median", "Max.")
  for (i in seq_along(nms)) {
    cat(sprintf("  %-8s  %8.*f  %8.*f  %8.*f\n",
                nms[i], digits, x$mean[i], digits, x$lower[i], digits, x$upper[i]))
  }

  cat("\n")
  invisible(x)
}

#' @export
print.ShrinkageTrees <- function(x, ...) {
  cat("\n")
  cat("ShrinkageTrees model\n")
  cat("---------------------\n")

  # Outcome type
  outcome_label <- switch(
    x$outcome_type,
    "continuous" = "Continuous",
    "binary" = "Binary (probit)",
    "right-censored" = paste0("Right-censored (AFT, timescale = ",
                              x$timescale, ")"),
    x$outcome_type
  )

  lbl <- function(s) sprintf("%-22s", s)

  cat(lbl("Outcome type:"), outcome_label, "\n", sep = "")
  cat(lbl("Prior:"), x$prior$prior_type_user, "\n", sep = "")
  cat(lbl("Number of trees:"), x$mcmc$number_of_trees, "\n", sep = "")
  cat(lbl("Training size (n):"), x$data_info$n_train, "\n", sep = "")
  cat(lbl("Number of features:"), x$data_info$p_features, "\n", sep = "")
  cat(lbl("Posterior draws:"), x$mcmc$N_post,
      " (burn-in ", x$mcmc$N_burn, ")\n", sep = "")

  if (!is.null(x$acceptance_ratio))
    cat(lbl("Acceptance ratio:"), round(mean(x$acceptance_ratio), 3), "\n",
        sep = "")

  if (!x$preprocess$sigma_known && !is.null(x$sigma))
    cat(lbl("Posterior mean sigma:"), round(mean(x$sigma), 3), "\n", sep = "")

  cat("\n")
  invisible(x)
}

# ── summary ──────────────────────────────────────────────────────────────────

.outcome_label <- function(outcome_type, timescale) {
  switch(
    outcome_type,
    "continuous"     = "Continuous",
    "binary"         = "Binary (probit)",
    "right-censored" = paste0("Right-censored (AFT, timescale = ",
                               timescale, ")"),
    outcome_type
  )
}

.sigma_summary <- function(sigma) {
  c(
    mean  = mean(sigma),
    sd    = sd(sigma),
    lower = unname(quantile(sigma, 0.025)),
    upper = unname(quantile(sigma, 0.975))
  )
}

.pred_summary <- function(preds) {
  c(mean = mean(preds), sd = sd(preds),
    min  = min(preds),  max = max(preds))
}

.vi_from_counts <- function(counts) {
  vi <- colMeans(counts / rowSums(counts))
  names(vi) <- if (!is.null(colnames(counts))) colnames(counts)
               else paste0("X", seq_along(vi))
  sort(vi, decreasing = TRUE)
}

.vi_from_probs <- function(probs) {
  vi <- colMeans(probs)
  names(vi) <- if (!is.null(colnames(probs))) colnames(probs)
               else paste0("X", seq_along(vi))
  sort(vi, decreasing = TRUE)
}

#' @export
summary.ShrinkageTrees <- function(object, ...) {

  out <- list(
    call         = object$call,
    outcome_type = object$outcome_type,
    timescale    = object$timescale,
    prior        = object$prior,
    mcmc         = object$mcmc,
    data_info    = object$data_info
  )

  if (!object$preprocess$sigma_known && !is.null(object$sigma))
    out$sigma <- .sigma_summary(object$sigma)

  out$predictions <- list(
    train = .pred_summary(object$train_predictions),
    test  = if (isTRUE(object$data_info$test_provided))
              .pred_summary(object$test_predictions)
  )

  if (isTRUE(object$prior$dirichlet) && !is.null(object$split_probs))
    out$variable_importance <- .vi_from_probs(object$split_probs)
  else if (!is.null(object$store_split_counts))
    out$variable_importance <- .vi_from_counts(object$store_split_counts)

  out$acceptance_ratio <- object$acceptance_ratio

  class(out) <- "summary.ShrinkageTrees"
  out
}

#' @export
print.summary.ShrinkageTrees <- function(x, n_vi = 10, ...) {

  cat("\n")
  cat("ShrinkageTrees model summary\n")
  cat("============================\n")
  cat("Call: "); print(x$call)
  cat("\n")
  cat("Outcome: ", .outcome_label(x$outcome_type, x$timescale),
      " | Prior: ", x$prior$prior_type_user,
      " | Trees: ", x$mcmc$number_of_trees, "\n", sep = "")
  cat("Data:    n = ", x$data_info$n_train,
      ", p = ", x$data_info$p_features,
      " | Draws: ", x$mcmc$N_post,
      " (burn-in ", x$mcmc$N_burn, ")\n", sep = "")

  if (!is.null(x$sigma)) {
    cat("\nPosterior sigma:\n")
    cat("  Mean: ", round(x$sigma["mean"],  3),
        "  SD: ",   round(x$sigma["sd"],    3),
        "  95% CI: [", round(x$sigma["lower"], 3),
        ", ",          round(x$sigma["upper"], 3), "]\n", sep = "")
  }

  cat("\nPredictions (posterior mean):\n")
  fmt_pred <- function(p) {
    paste0("mean = ", round(p["mean"], 3),
           ", sd = ",  round(p["sd"],   3),
           ", range = [", round(p["min"], 3), ", ", round(p["max"], 3), "]")
  }
  cat("  Train: ", fmt_pred(x$predictions$train), "\n", sep = "")
  if (!is.null(x$predictions$test))
    cat("  Test:  ", fmt_pred(x$predictions$test), "\n", sep = "")

  if (!is.null(x$variable_importance)) {
    cat("\nVariable importance (posterior inclusion probability):\n")
    vi <- head(x$variable_importance, n_vi)
    cat(" ", paste(names(vi), round(vi, 3), sep = ": ", collapse = "   "), "\n")
  }

  if (!is.null(x$acceptance_ratio))
    cat("\nMCMC acceptance ratio: ", round(mean(x$acceptance_ratio), 3), "\n",
        sep = "")

  cat("\n")
  invisible(x)
}

#' @export
summary.CausalShrinkageForest <- function(object, ...) {


  out <- list(
    call         = object$call,
    outcome_type = object$outcome_type,
    timescale    = object$timescale,
    prior        = object$prior,
    mcmc         = object$mcmc,
    data_info    = object$data_info
  )

  tau <- object$train_predictions_treat
  te  <- list(ate = mean(tau), cate_sd = sd(tau))

  if (!is.null(object$train_predictions_sample_treat)) {
    ate_samples  <- rowMeans(object$train_predictions_sample_treat)
    te$ate_lower <- unname(quantile(ate_samples, 0.025))
    te$ate_upper <- unname(quantile(ate_samples, 0.975))
  }
  out$treatment_effect <- te

  out$prognostic <- .pred_summary(object$train_predictions_control)

  if (!object$preprocess$sigma_known && !is.null(object$sigma))
    out$sigma <- .sigma_summary(object$sigma)

  if (isTRUE(object$prior$control$dirichlet) && !is.null(object$split_probs_control))
    out$variable_importance_control <- .vi_from_probs(object$split_probs_control)
  else if (!is.null(object$split_counts_control))
    out$variable_importance_control <- .vi_from_counts(object$split_counts_control)

  if (isTRUE(object$prior$treat$dirichlet) && !is.null(object$split_probs_treat))
    out$variable_importance_treat <- .vi_from_probs(object$split_probs_treat)
  else if (!is.null(object$split_counts_treat))
    out$variable_importance_treat <- .vi_from_counts(object$split_counts_treat)

  out$acceptance_ratios <- list(
    control = object$acceptance_ratio_control,
    treat   = object$acceptance_ratio_treat
  )

  class(out) <- "summary.CausalShrinkageForest"
  out
}

#' @export
print.summary.CausalShrinkageForest <- function(x, n_vi = 10, ...) {

  cat("\n")
  cat("CausalShrinkageForest model summary\n")
  cat("=====================================\n")
  cat("Call: "); print(x$call)
  cat("\n")
  cat("Outcome: ", .outcome_label(x$outcome_type, x$timescale), "\n", sep = "")
  cat("Prior:   control = ", x$prior$control$prior_type_user,
      ", treatment = ", x$prior$treat$prior_type_user, "\n", sep = "")
  cat("Trees:   control = ", x$mcmc$number_of_trees_control,
      ", treatment = ", x$mcmc$number_of_trees_treat, "\n", sep = "")
  cat("Data:    n = ", x$data_info$n_train,
      ", p_control = ", x$data_info$p_control,
      ", p_treat = ", x$data_info$p_treat,
      " | Draws: ", x$mcmc$N_post,
      " (burn-in ", x$mcmc$N_burn, ")\n", sep = "")

  cat("\nTreatment effect (CATE):\n")
  te <- x$treatment_effect
  if (!is.null(te$ate_lower)) {
    cat("  ATE:     ", round(te$ate, 4),
        "  95% CI: [", round(te$ate_lower, 4), ", ",
                       round(te$ate_upper, 4), "]\n", sep = "")
  } else {
    cat("  ATE:     ", round(te$ate, 4),
        "  (no CI — refit with store_posterior_sample = TRUE)\n", sep = "")
  }
  cat("  CATE SD: ", round(te$cate_sd, 4), "\n", sep = "")

  cat("\nPrognostic function (mu):\n")
  p <- x$prognostic
  cat("  Mean: ", round(p["mean"], 3),
      "  SD: ",   round(p["sd"],   3),
      "  Range: [", round(p["min"], 3), ", ", round(p["max"], 3), "]\n",
      sep = "")

  if (!is.null(x$sigma)) {
    cat("\nPosterior sigma:\n")
    cat("  Mean: ", round(x$sigma["mean"],  3),
        "  SD: ",   round(x$sigma["sd"],    3),
        "  95% CI: [", round(x$sigma["lower"], 3),
        ", ",          round(x$sigma["upper"], 3), "]\n", sep = "")
  }

  fmt_vi <- function(vi) {
    vi <- head(vi, n_vi)
    paste(names(vi), round(vi, 3), sep = ": ", collapse = "   ")
  }
  if (!is.null(x$variable_importance_control)) {
    cat("\nVariable importance — control forest (posterior inclusion probability):\n")
    cat(" ", fmt_vi(x$variable_importance_control), "\n")
  }
  if (!is.null(x$variable_importance_treat)) {
    cat("\nVariable importance — treatment forest (posterior inclusion probability):\n")
    cat(" ", fmt_vi(x$variable_importance_treat), "\n")
  }

  cat("\nMCMC acceptance ratios:",
      " control = ", round(x$acceptance_ratios$control, 3),
      ", treatment = ", round(x$acceptance_ratios$treat, 3), "\n", sep = "")

  cat("\n")
  invisible(x)
}

#' @export
print.CausalShrinkageForest <- function(x, ...) {

  cat("\n")
  cat("CausalShrinkageForest model\n")
  cat("---------------------------\n")

  lbl <- function(s) sprintf("%-22s", s)

  # ── General info ──────────────────────────────────────────────────────────
  cat(lbl("Outcome type:"), .outcome_label(x$outcome_type, x$timescale),
      "\n", sep = "")
  cat(lbl("Training size (n):"), x$data_info$n_train, "\n", sep = "")
  cat(lbl("Posterior draws:"), x$mcmc$N_post,
      " (burn-in ", x$mcmc$N_burn, ")\n", sep = "")

  if (!x$preprocess$sigma_known && !is.null(x$sigma))
    cat(lbl("Posterior mean sigma:"), round(mean(x$sigma), 3), "\n", sep = "")

  # ── Per-forest info (side by side) ────────────────────────────────────────
  col  <- function(s) sprintf("%-20s", s)
  hdr  <- function(s) sprintf("%-20s", s)

  cat("\n")
  cat(sprintf("%-22s", ""), hdr("Control"), hdr("Treatment"), "\n", sep = "")
  cat(sprintf("%-22s", ""), strrep("-", 19), " ", strrep("-", 19), "\n",
      sep = "")

  cat(lbl("Prior:"),
      col(x$prior$control$prior_type_user),
      col(x$prior$treat$prior_type_user), "\n", sep = "")
  cat(lbl("Number of trees:"),
      col(x$mcmc$number_of_trees_control),
      col(x$mcmc$number_of_trees_treat), "\n", sep = "")
  cat(lbl("Number of features:"),
      col(x$data_info$p_control),
      col(x$data_info$p_treat), "\n", sep = "")

  if (!is.null(x$acceptance_ratio_control) ||
      !is.null(x$acceptance_ratio_treat)) {
    acc_c <- if (!is.null(x$acceptance_ratio_control))
               round(mean(x$acceptance_ratio_control), 3) else "-"
    acc_t <- if (!is.null(x$acceptance_ratio_treat))
               round(mean(x$acceptance_ratio_treat), 3) else "-"
    cat(lbl("Acceptance ratio:"), col(acc_c), col(acc_t), "\n", sep = "")
  }

  cat("\n")
  invisible(x)
}