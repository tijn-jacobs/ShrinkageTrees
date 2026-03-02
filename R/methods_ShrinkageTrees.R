#' @export
print.ShrinkageTrees <- function(x, ...) {

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

  cat("Outcome type:      ", outcome_label, "\n", sep = "")
  cat("Prior:             ", x$prior$prior_type_user, "\n", sep = "")
  cat("Number of trees:   ", x$mcmc$number_of_trees, "\n", sep = "")
  cat("Training size (n): ", x$data_info$n_train, "\n", sep = "")
  cat("Number of features:", x$data_info$p_features, "\n", sep = "")
  cat("Posterior draws:   ", x$mcmc$N_post,
      " (burn-in ", x$mcmc$N_burn, ")\n", sep = "")

  if (!is.null(x$acceptance_ratio)) {
    cat("Acceptance ratio:  ",
        round(mean(x$acceptance_ratio), 3), "\n", sep = "")
  }

  if (!x$preprocess$sigma_known && !is.null(x$sigma)) {
    cat("Posterior mean sigma: ",
        round(mean(x$sigma), 3), "\n", sep = "")
  }

  invisible(x)
}

#' @export
print.CausalShrinkageForest <- function(x, ...) {

  cat("CausalShrinkageForest model\n")
  cat("---------------------------\n")

  outcome_label <- switch(
    x$outcome_type,
    "continuous" = "Continuous",
    "binary" = "Binary (probit)",
    "right-censored" = paste0("Right-censored (AFT, timescale = ",
                              x$timescale, ")"),
    x$outcome_type
  )

  cat("Outcome type:         ", outcome_label, "\n", sep = "")
  cat("Prior (control):      ", x$prior$control$prior_type_user, "\n", sep = "")
  cat("Prior (treatment):    ", x$prior$treat$prior_type_user, "\n", sep = "")
  cat("Trees (control):      ", x$mcmc$number_of_trees_control, "\n", sep = "")
  cat("Trees (treatment):    ", x$mcmc$number_of_trees_treat, "\n", sep = "")
  cat("Training size (n):    ", x$data_info$n_train, "\n", sep = "")
  cat("Features (control):   ", x$data_info$p_control, "\n", sep = "")
  cat("Features (treatment): ", x$data_info$p_treat, "\n", sep = "")
  cat("Posterior draws:      ", x$mcmc$N_post,
      " (burn-in ", x$mcmc$N_burn, ")\n", sep = "")

  if (!is.null(x$acceptance_ratio_control)) {
    cat("Acceptance (control): ",
        round(mean(x$acceptance_ratio_control), 3), "\n", sep = "")
  }

  if (!is.null(x$acceptance_ratio_treat)) {
    cat("Acceptance (treatment):",
        round(mean(x$acceptance_ratio_treat), 3), "\n", sep = "")
  }

  if (!x$preprocess$sigma_known && !is.null(x$sigma)) {
    cat("Posterior mean sigma: ",
        round(mean(x$sigma), 3), "\n", sep = "")
  }

  invisible(x)
}