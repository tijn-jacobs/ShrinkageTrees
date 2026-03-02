#' Posterior predictive inference for a ShrinkageTrees model
#'
#' Re-runs the MCMC sampler on new covariate data using the stored training
#' data and hyperparameters, returning posterior mean predictions and credible
#' interval bounds.
#'
#' @param object A fitted \code{ShrinkageTrees} model object.
#' @param newdata A matrix (or object coercible to one) of new covariates with
#'   the same number of columns as the training data.
#' @param level Credible interval width. Default \code{0.95}.
#' @param ... Currently unused.
#' @return A \code{ShrinkageTreesPrediction} object with elements:
#'   \describe{
#'     \item{mean}{Posterior mean predictions (length \code{nrow(newdata)}).}
#'     \item{lower}{Lower credible interval bound.}
#'     \item{upper}{Upper credible interval bound.}
#'     \item{n}{Number of test observations.}
#'     \item{level}{Credible level used.}
#'     \item{outcome_type}{Outcome type inherited from the fitted model.}
#'     \item{timescale}{Timescale inherited from the fitted model (survival only).}
#'   }
#' @seealso \code{\link{HorseTrees}}, \code{\link{ShrinkageTrees}},
#'   \code{\link{print.ShrinkageTreesPrediction}},
#'   \code{\link{summary.ShrinkageTreesPrediction}}
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

#' Print a ShrinkageTreesPrediction object
#'
#' Displays a formatted table of posterior mean predictions and credible
#' interval bounds for the first \code{n_head} observations.
#'
#' @param x A \code{ShrinkageTreesPrediction} object.
#' @param n_head Number of observations to display. Default \code{6}.
#' @param digits Number of decimal places. Default \code{3}.
#' @param ... Currently unused.
#' @return Invisibly returns \code{x}.
#' @seealso \code{\link{predict.ShrinkageTrees}},
#'   \code{\link{summary.ShrinkageTreesPrediction}}
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

#' Summarise a ShrinkageTreesPrediction object
#'
#' Returns distributional summaries (min, Q1, median, max) of the posterior
#' mean predictions and credible interval bounds across all observations.
#'
#' @param object A \code{ShrinkageTreesPrediction} object.
#' @param ... Currently unused.
#' @return A \code{summary.ShrinkageTreesPrediction} object.
#' @seealso \code{\link{predict.ShrinkageTrees}},
#'   \code{\link{print.summary.ShrinkageTreesPrediction}}
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

#' Print a ShrinkageTreesPrediction summary
#'
#' Displays distributional summaries (min, Q1, median, max) of the posterior
#' mean predictions and credible interval bounds.
#'
#' @param x A \code{summary.ShrinkageTreesPrediction} object.
#' @param digits Number of decimal places. Default \code{3}.
#' @param ... Currently unused.
#' @return Invisibly returns \code{x}.
#' @seealso \code{\link{summary.ShrinkageTreesPrediction}}
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

#' Print a ShrinkageTrees model
#'
#' Displays a concise summary of a fitted \code{ShrinkageTrees} model,
#' including outcome type, prior, MCMC settings, acceptance ratio, and
#' posterior mean sigma.
#'
#' @param x A fitted \code{ShrinkageTrees} model object.
#' @param ... Currently unused.
#' @return Invisibly returns \code{x}.
#' @seealso \code{\link{summary.ShrinkageTrees}}, \code{\link{HorseTrees}},
#'   \code{\link{ShrinkageTrees}}
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
  n_ch <- x$mcmc$n_chains
  if (!is.null(n_ch) && n_ch > 1) {
    cat(lbl("Chains:"),          n_ch,         "\n", sep = "")
    cat(lbl("Draws per chain:"), x$mcmc$N_post,
        " (burn-in ", x$mcmc$N_burn, ")\n", sep = "")
  } else {
    cat(lbl("Posterior draws:"), x$mcmc$N_post,
        " (burn-in ", x$mcmc$N_burn, ")\n", sep = "")
  }

  if (!is.null(x$chains)) {
    cat(lbl("Acceptance ratio:"),
        paste(round(x$chains$acceptance_ratios, 3), collapse = ", "),
        " (per chain)\n", sep = "")
  } else if (!is.null(x$acceptance_ratio)) {
    cat(lbl("Acceptance ratio:"), round(mean(x$acceptance_ratio), 3), "\n",
        sep = "")
  }

  if (!x$preprocess$sigma_known && !is.null(x$sigma))
    cat(lbl("Posterior mean sigma:"), round(mean(x$sigma), 3), "\n", sep = "")

  cat("\n")
  invisible(x)
}

# -- summary ------------------------------------------------------------------

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

#' Summarise a ShrinkageTrees model
#'
#' Returns an inspectable list with posterior sigma summaries, prediction
#' summaries, variable importance (posterior inclusion probabilities), and
#' MCMC diagnostics.
#'
#' @param object A fitted \code{ShrinkageTrees} model object.
#' @param ... Currently unused.
#' @return A \code{summary.ShrinkageTrees} object with elements:
#'   \describe{
#'     \item{call}{The original model call.}
#'     \item{outcome_type}{Outcome type (\code{"continuous"}, \code{"binary"},
#'       or \code{"right-censored"}).}
#'     \item{timescale}{Timescale for survival outcomes (\code{"time"} or
#'       \code{"log"}).}
#'     \item{prior}{Prior specification.}
#'     \item{mcmc}{MCMC settings (trees, draws, burn-in).}
#'     \item{data_info}{Training and test data dimensions.}
#'     \item{sigma}{Named vector with posterior mean, SD, and 95\% CI of sigma
#'       (continuous and survival outcomes only).}
#'     \item{predictions}{List with \code{train} (and optionally \code{test})
#'       prediction summaries (mean, SD, range).}
#'     \item{variable_importance}{Named vector of posterior inclusion
#'       probabilities, sorted decreasingly (if available).}
#'     \item{acceptance_ratio}{MCMC acceptance ratio vector.}
#'   }
#' @seealso \code{\link{print.summary.ShrinkageTrees}},
#'   \code{\link{HorseTrees}}, \code{\link{ShrinkageTrees}}
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
  if (!is.null(object$chains))
    out$chains <- object$chains

  class(out) <- "summary.ShrinkageTrees"
  out
}

#' Print a ShrinkageTrees model summary
#'
#' Displays a detailed summary of a \code{ShrinkageTrees} model, including
#' model specification, posterior sigma, prediction summaries, variable
#' importance, and MCMC diagnostics.
#'
#' @param x A \code{summary.ShrinkageTrees} object.
#' @param n_vi Maximum number of variables to display in the variable
#'   importance table. Default \code{10}.
#' @param ... Currently unused.
#' @return Invisibly returns \code{x}.
#' @seealso \code{\link{summary.ShrinkageTrees}}
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
  n_ch <- x$mcmc$n_chains
  draws_str <- if (!is.null(n_ch) && n_ch > 1)
    paste0(x$mcmc$N_post, " x ", n_ch, " chains (burn-in ", x$mcmc$N_burn, ")")
  else
    paste0(x$mcmc$N_post, " (burn-in ", x$mcmc$N_burn, ")")
  cat("Data:    n = ", x$data_info$n_train,
      ", p = ", x$data_info$p_features,
      " | Draws: ", draws_str, "\n", sep = "")

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

  if (!is.null(x$chains)) {
    cat("\nMCMC acceptance ratio (per chain): ",
        paste(round(x$chains$acceptance_ratios, 3), collapse = ", "), "\n",
        sep = "")
  } else if (!is.null(x$acceptance_ratio)) {
    cat("\nMCMC acceptance ratio: ", round(mean(x$acceptance_ratio), 3), "\n",
        sep = "")
  }

  cat("\n")
  invisible(x)
}

#' Summarise a CausalShrinkageForest model
#'
#' Returns an inspectable list with treatment effect estimates, prognostic
#' function summaries, posterior sigma, variable importance for each forest,
#' and MCMC diagnostics.
#'
#' @param object A fitted \code{CausalShrinkageForest} model object.
#' @param ... Currently unused.
#' @return A \code{summary.CausalShrinkageForest} object with elements:
#'   \describe{
#'     \item{call}{The original model call.}
#'     \item{outcome_type}{Outcome type.}
#'     \item{timescale}{Timescale for survival outcomes.}
#'     \item{prior}{Prior specification for control and treatment forests.}
#'     \item{mcmc}{MCMC settings.}
#'     \item{data_info}{Training and test data dimensions.}
#'     \item{treatment_effect}{List with \code{ate} (posterior mean ATE),
#'       \code{cate_sd} (SD of individual CATEs), and optionally
#'       \code{ate_lower} and \code{ate_upper} (95\% CI; requires
#'       \code{store_posterior_sample = TRUE}).}
#'     \item{prognostic}{Summary of the prognostic function (mean, SD, range).}
#'     \item{sigma}{Named vector with posterior mean, SD, and 95\% CI of sigma
#'       (if estimated).}
#'     \item{variable_importance_control}{Variable importance for the control
#'       forest (if available).}
#'     \item{variable_importance_treat}{Variable importance for the treatment
#'       forest (if available).}
#'     \item{acceptance_ratios}{List with acceptance ratios for each forest.}
#'   }
#' @seealso \code{\link{print.summary.CausalShrinkageForest}},
#'   \code{\link{CausalHorseForest}}, \code{\link{CausalShrinkageForest}}
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
  if (!is.null(object$chains))
    out$acceptance_ratios$per_chain <- object$chains

  class(out) <- "summary.CausalShrinkageForest"
  out
}

#' Print a CausalShrinkageForest model summary
#'
#' Displays a detailed summary of a \code{CausalShrinkageForest} model,
#' including model specification, treatment effect estimates, prognostic
#' function, posterior sigma, variable importance for each forest, and MCMC
#' diagnostics.
#'
#' @param x A \code{summary.CausalShrinkageForest} object.
#' @param n_vi Maximum number of variables to display per variable importance
#'   table. Default \code{10}.
#' @param ... Currently unused.
#' @return Invisibly returns \code{x}.
#' @seealso \code{\link{summary.CausalShrinkageForest}}
#' @export
print.summary.CausalShrinkageForest <- function(x, n_vi = 10, ...) {
#TESTESTtest
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
  n_ch <- x$mcmc$n_chains
  draws_str <- if (!is.null(n_ch) && n_ch > 1)
    paste0(x$mcmc$N_post, " x ", n_ch, " chains (burn-in ", x$mcmc$N_burn, ")")
  else
    paste0(x$mcmc$N_post, " (burn-in ", x$mcmc$N_burn, ")")
  cat("Data:    n = ", x$data_info$n_train,
      ", p_control = ", x$data_info$p_control,
      ", p_treat = ", x$data_info$p_treat,
      " | Draws: ", draws_str, "\n", sep = "")

  cat("\nTreatment effect (CATE):\n")
  te <- x$treatment_effect
  if (!is.null(te$ate_lower)) {
    cat("  ATE:     ", round(te$ate, 4),
        "  95% CI: [", round(te$ate_lower, 4), ", ",
                       round(te$ate_upper, 4), "]\n", sep = "")
  } else {
    cat("  ATE:     ", round(te$ate, 4),
        "  (no CI -- refit with store_posterior_sample = TRUE)\n", sep = "")
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
    cat("\nVariable importance - control forest (posterior inclusion probability):\n")
    cat(" ", fmt_vi(x$variable_importance_control), "\n")
  }
  if (!is.null(x$variable_importance_treat)) {
    cat("\nVariable importance - treatment forest (posterior inclusion probability):\n")
    cat(" ", fmt_vi(x$variable_importance_treat), "\n")
  }

  if (!is.null(x$acceptance_ratios$per_chain)) {
    pc <- x$acceptance_ratios$per_chain
    cat("\nMCMC acceptance ratios (per chain):\n")
    cat("  control:   ",
        paste(round(pc$acceptance_ratios_control, 3), collapse = ", "), "\n",
        sep = "")
    cat("  treatment: ",
        paste(round(pc$acceptance_ratios_treat,   3), collapse = ", "), "\n",
        sep = "")
  } else {
    cat("\nMCMC acceptance ratios:",
        " control = ", round(mean(x$acceptance_ratios$control), 3),
        ", treatment = ", round(mean(x$acceptance_ratios$treat), 3), "\n",
        sep = "")
  }

  cat("\n")
  invisible(x)
}

#' Print a CausalShrinkageForest model
#'
#' Displays a concise summary of a fitted \code{CausalShrinkageForest} model
#' with per-forest columns for priors, tree counts, feature counts, and MCMC
#' acceptance ratios.
#'
#' @param x A fitted \code{CausalShrinkageForest} model object.
#' @param ... Currently unused.
#' @return Invisibly returns \code{x}.
#' @seealso \code{\link{summary.CausalShrinkageForest}},
#'   \code{\link{CausalHorseForest}}, \code{\link{CausalShrinkageForest}}
#' @export
print.CausalShrinkageForest <- function(x, ...) {

  cat("\n")
  cat("CausalShrinkageForest model\n")
  cat("---------------------------\n")

  lbl <- function(s) sprintf("%-22s", s)

  # -- General info ----------------------------------------------------------
  cat(lbl("Outcome type:"), .outcome_label(x$outcome_type, x$timescale),
      "\n", sep = "")
  cat(lbl("Training size (n):"), x$data_info$n_train, "\n", sep = "")
  n_ch <- x$mcmc$n_chains
  if (!is.null(n_ch) && n_ch > 1) {
    cat(lbl("Chains:"),          n_ch,          "\n", sep = "")
    cat(lbl("Draws per chain:"), x$mcmc$N_post,
        " (burn-in ", x$mcmc$N_burn, ")\n", sep = "")
  } else {
    cat(lbl("Posterior draws:"), x$mcmc$N_post,
        " (burn-in ", x$mcmc$N_burn, ")\n", sep = "")
  }

  if (!x$preprocess$sigma_known && !is.null(x$sigma))
    cat(lbl("Posterior mean sigma:"), round(mean(x$sigma), 3), "\n", sep = "")

  # -- Per-forest info (side by side) ----------------------------------------
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

  if (!is.null(x$chains)) {
    acc_c <- paste(round(x$chains$acceptance_ratios_control, 3), collapse = ", ")
    acc_t <- paste(round(x$chains$acceptance_ratios_treat,   3), collapse = ", ")
    cat(lbl("Acceptance ratio:"), col(acc_c), col(acc_t), "(per chain)\n",
        sep = "")
  } else if (!is.null(x$acceptance_ratio_control) ||
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