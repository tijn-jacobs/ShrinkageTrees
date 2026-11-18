# Scenarios for comparing the refactored code against the code before it.
#
# A baseline was recorded by running every scenario with the OLD package
# (commit 8d4bc0d) and saving the results to
# tests/testthat/fixtures/refactor-baseline.rds. That file is kept locally and
# not committed. tests/testthat/test-refactor-equivalence.R runs the scenarios
# with the CURRENT package and compares.
#
# Rules for a scenario:
#   * It sets its own seed, builds its own data and returns a fitted object.
#   * It uses exported functions only, so it runs against either version.
#   * `baseline = TRUE` means "running on the old code". Only the two
#     horseshoe_fw scenarios use it: the old code is called with
#     prior_type = "horseshoe_fw", the new code with "half-cauchy".

refactor_data <- function(seed, n = 40, p = 4) {
  set.seed(seed)
  X      <- matrix(runif(n * p), n, p)
  X_test <- matrix(runif(10 * p), 10, p)
  trt    <- rbinom(n, 1, 0.5)
  mu     <- 2 * X[, 1] - X[, 2]

  y_cont <- as.numeric(mu + rnorm(n, 0, 0.5))
  y_bin  <- as.numeric(mu + rnorm(n) > 0.5)

  # Right-censored survival times
  t_event <- exp(0.5 * mu + trt * 0.4 + rnorm(n, 0, 0.4))
  t_cens  <- exp(rnorm(n, 0.8, 0.5))
  time    <- pmin(t_event, t_cens)
  status  <- as.numeric(t_event <= t_cens)

  # Interval-censored version: a third exact, a third right-censored,
  # a third observed in a finite interval
  kind  <- rep(1:3, length.out = n)
  left  <- ifelse(kind == 3, t_event * 0.8, t_event)
  right <- ifelse(kind == 1, t_event, ifelse(kind == 2, Inf, t_event * 1.3))

  list(X = X, X_test = X_test, trt = trt, y_cont = y_cont, y_bin = y_bin,
       y_causal = as.numeric(y_cont + trt * (1 + X[, 3])),
       time = time, status = status, left = left, right = right,
       propensity = rep(0.5, n), n = n, p = p)
}

refactor_scenarios <- function() {
  mcmc <- list(N_post = 20, N_burn = 10, verbose = FALSE)
  fit  <- function(fun, ...) do.call(fun, c(list(...), mcmc))

  list(

    # ── Single forest: HorseTrees, one per outcome type ──────────────────────
    horsetrees_continuous = function(baseline = FALSE) {
      d <- refactor_data(101)
      fit(HorseTrees, y = d$y_cont, X_train = d$X, X_test = d$X_test,
          outcome_type = "continuous", number_of_trees = 5)
    },
    horsetrees_binary = function(baseline = FALSE) {
      d <- refactor_data(102)
      fit(HorseTrees, y = d$y_bin, X_train = d$X, X_test = d$X_test,
          outcome_type = "binary", number_of_trees = 5)
    },
    horsetrees_right_censored = function(baseline = FALSE) {
      d <- refactor_data(103)
      fit(HorseTrees, y = d$time, status = d$status, X_train = d$X,
          X_test = d$X_test, outcome_type = "right-censored",
          timescale = "time", number_of_trees = 5, k = 0.8)
    },
    horsetrees_interval_censored = function(baseline = FALSE) {
      d <- refactor_data(104)
      fit(HorseTrees, left_time = d$left, right_time = d$right,
          X_train = d$X, outcome_type = "interval-censored",
          timescale = "time", number_of_trees = 5)
    },
    horsetrees_no_samples = function(baseline = FALSE) {
      d <- refactor_data(105)
      fit(HorseTrees, y = d$y_cont, X_train = d$X,
          outcome_type = "continuous", number_of_trees = 5,
          store_posterior_sample = FALSE)
    },

    # ── Single forest: ShrinkageTrees, one per prior ─────────────────────────
    shrinkagetrees_horseshoe_default = function(baseline = FALSE) {
      d <- refactor_data(111)
      fit(ShrinkageTrees, y = d$y_cont, X_train = d$X, X_test = d$X_test,
          outcome_type = "continuous", prior_type = "horseshoe",
          number_of_trees = 5)
    },
    shrinkagetrees_halfcauchy = function(baseline = FALSE) {
      d <- refactor_data(112)
      fit(ShrinkageTrees, y = d$y_cont, X_train = d$X, X_test = d$X_test,
          outcome_type = "continuous", prior_type = "half-cauchy",
          local_hp = 0.3, number_of_trees = 5)
    },
    shrinkagetrees_standard_log_time = function(baseline = FALSE) {
      d <- refactor_data(113)
      fit(ShrinkageTrees, y = log(d$time), status = d$status, X_train = d$X,
          outcome_type = "right-censored", timescale = "log",
          prior_type = "standard", local_hp = 0.2, number_of_trees = 5)
    },
    shrinkagetrees_dirichlet = function(baseline = FALSE) {
      d <- refactor_data(114)
      fit(ShrinkageTrees, y = d$y_cont, X_train = d$X, X_test = d$X_test,
          outcome_type = "continuous", prior_type = "dirichlet",
          local_hp = 0.2, number_of_trees = 5)
    },
    shrinkagetrees_known_sigma = function(baseline = FALSE) {
      d <- refactor_data(115)
      fit(ShrinkageTrees, y = d$y_cont, X_train = d$X,
          outcome_type = "continuous", prior_type = "horseshoe",
          sigma = 0.7, number_of_trees = 5)
    },
    shrinkagetrees_interval_log = function(baseline = FALSE) {
      d <- refactor_data(116)
      fit(ShrinkageTrees, left_time = log(d$left), right_time = log(d$right),
          X_train = d$X, outcome_type = "interval-censored",
          timescale = "log", prior_type = "horseshoe", number_of_trees = 5)
    },
    shrinkagetrees_two_chains = function(baseline = FALSE) {
      # One core makes mclapply() run the chains in sequence, so the seed
      # fixes the result. This exercises the multi-chain argument passing.
      old <- Sys.getenv("OMP_THREAD_LIMIT", unset = NA)
      Sys.setenv(OMP_THREAD_LIMIT = "1")
      on.exit(if (is.na(old)) Sys.unsetenv("OMP_THREAD_LIMIT")
              else Sys.setenv(OMP_THREAD_LIMIT = old))
      d <- refactor_data(117)
      fit(ShrinkageTrees, y = d$y_cont, X_train = d$X, X_test = d$X_test,
          outcome_type = "continuous", prior_type = "horseshoe",
          number_of_trees = 5, n_chains = 2)
    },

    # ── Removed option: old "horseshoe_fw" against new "half-cauchy" ─────────
    # The stored hyperparameters are dropped before comparing: the old code
    # recorded the global scale it was given (and then ignored), the new code
    # records the placeholder 1.
    horseshoe_fw_single = function(baseline = FALSE) {
      d <- refactor_data(121)
      f <- if (baseline) {
        fit(ShrinkageTrees, y = d$y_cont, X_train = d$X, X_test = d$X_test,
            outcome_type = "continuous", prior_type = "horseshoe_fw",
            local_hp = 0.3, global_hp = 0.3, number_of_trees = 5)
      } else {
        fit(ShrinkageTrees, y = d$y_cont, X_train = d$X, X_test = d$X_test,
            outcome_type = "continuous", prior_type = "half-cauchy",
            local_hp = 0.3, number_of_trees = 5)
      }
      f$args <- NULL
      f
    },
    horseshoe_fw_causal = function(baseline = FALSE) {
      d <- refactor_data(122)
      prior <- if (baseline) "horseshoe_fw" else "half-cauchy"
      hp    <- if (baseline) 0.3 else NULL
      f <- fit(CausalShrinkageForest, y = d$y_causal, X_train_control = d$X,
               X_train_treat = d$X, treatment_indicator_train = d$trt,
               outcome_type = "continuous",
               prior_type_control = prior, prior_type_treat = prior,
               local_hp_control = 0.3, local_hp_treat = 0.3,
               global_hp_control = hp, global_hp_treat = hp,
               number_of_trees_control = 5, number_of_trees_treat = 5,
               store_posterior_sample = TRUE)
      f$args <- NULL
      f
    },

    # ── Survival wrappers ────────────────────────────────────────────────────
    survivalbart_right_censored = function(baseline = FALSE) {
      d <- refactor_data(131)
      fit(SurvivalBART, time = d$time, status = d$status, X_train = d$X,
          X_test = d$X_test, number_of_trees = 5)
    },
    survivaldart_interval_censored = function(baseline = FALSE) {
      d <- refactor_data(132)
      fit(SurvivalDART, left_time = d$left, right_time = d$right,
          X_train = d$X, number_of_trees = 5)
    },
    survivalbcf_right_censored = function(baseline = FALSE) {
      d <- refactor_data(133)
      fit(SurvivalBCF, time = d$time, status = d$status, X_train = d$X,
          treatment = d$trt, propensity = d$propensity,
          number_of_trees_control = 5, number_of_trees_treat = 5)
    },
    survivalshrinkagebcf_interval_censored = function(baseline = FALSE) {
      d <- refactor_data(134)
      fit(SurvivalShrinkageBCF, left_time = d$left, right_time = d$right,
          X_train = d$X, treatment = d$trt,
          number_of_trees_control = 5, number_of_trees_treat = 5)
    },

    # ── Causal: CausalHorseForest, one per outcome type ──────────────────────
    causalhorseforest_continuous = function(baseline = FALSE) {
      d <- refactor_data(141)
      fit(CausalHorseForest, y = d$y_causal, X_train_control = d$X,
          X_train_treat = d$X, treatment_indicator_train = d$trt,
          X_test_control = d$X_test, X_test_treat = d$X_test,
          outcome_type = "continuous", number_of_trees = 5,
          store_posterior_sample = TRUE)
    },
    causalhorseforest_right_censored = function(baseline = FALSE) {
      d <- refactor_data(142)
      fit(CausalHorseForest, y = d$time, status = d$status,
          X_train_control = d$X, X_train_treat = d$X,
          treatment_indicator_train = d$trt,
          outcome_type = "right-censored", timescale = "time",
          number_of_trees = 5, store_posterior_sample = TRUE)
    },
    causalhorseforest_interval_censored = function(baseline = FALSE) {
      d <- refactor_data(143)
      fit(CausalHorseForest, left_time = d$left, right_time = d$right,
          X_train_control = d$X, X_train_treat = d$X,
          treatment_indicator_train = d$trt,
          outcome_type = "interval-censored", timescale = "time",
          number_of_trees = 5, k = 1.2)
    },

    # ── Causal: CausalShrinkageForest, priors and treatment codings ──────────
    causalshrinkage_standard_binary_coding = function(baseline = FALSE) {
      d <- refactor_data(151)
      fit(CausalShrinkageForest, y = d$y_causal, X_train_control = d$X,
          X_train_treat = d$X, treatment_indicator_train = d$trt,
          outcome_type = "continuous",
          prior_type_control = "standard", prior_type_treat = "standard",
          local_hp_control = 0.2, local_hp_treat = 0.2,
          number_of_trees_control = 5, number_of_trees_treat = 5,
          treatment_coding = "binary", store_posterior_sample = TRUE)
    },
    causalshrinkage_halfcauchy_scales = function(baseline = FALSE) {
      d <- refactor_data(152)
      fit(CausalShrinkageForest, y = d$y_causal, X_train_control = d$X,
          X_train_treat = d$X, treatment_indicator_train = d$trt,
          outcome_type = "continuous",
          prior_type_control = "standard-halfcauchy",
          prior_type_treat = "dirichlet",
          local_hp_control = 1.5, local_hp_treat = 0.2,
          number_of_trees_control = 5, number_of_trees_treat = 5)
    },
    causalshrinkage_invariant_coding = function(baseline = FALSE) {
      d <- refactor_data(153)
      fit(CausalShrinkageForest, y = d$y_causal, X_train_control = d$X,
          X_train_treat = d$X, treatment_indicator_train = d$trt,
          outcome_type = "continuous",
          number_of_trees_control = 5, number_of_trees_treat = 3,
          treatment_coding = "invariant", store_posterior_sample = TRUE)
    },
    causalshrinkage_adaptive_coding_known_sigma = function(baseline = FALSE) {
      d <- refactor_data(154)
      fit(CausalShrinkageForest, y = d$y_causal, X_train_control = d$X,
          X_train_treat = d$X, treatment_indicator_train = d$trt,
          outcome_type = "continuous", sigma = 0.8,
          number_of_trees_control = 5, number_of_trees_treat = 5,
          treatment_coding = "adaptive", propensity = d$propensity)
    },
    causalhorseforest_two_chains = function(baseline = FALSE) {
      old <- Sys.getenv("OMP_THREAD_LIMIT", unset = NA)
      Sys.setenv(OMP_THREAD_LIMIT = "1")
      on.exit(if (is.na(old)) Sys.unsetenv("OMP_THREAD_LIMIT")
              else Sys.setenv(OMP_THREAD_LIMIT = old))
      d <- refactor_data(155)
      fit(CausalHorseForest, y = d$y_causal, X_train_control = d$X,
          X_train_treat = d$X, treatment_indicator_train = d$trt,
          outcome_type = "continuous", number_of_trees = 5,
          store_posterior_sample = TRUE, n_chains = 2)
    },

    # ── predict(): refits from the stored data and hyperparameters ───────────
    predict_horsetrees_continuous = function(baseline = FALSE) {
      d <- refactor_data(161)
      f <- fit(HorseTrees, y = d$y_cont, X_train = d$X,
               outcome_type = "continuous", number_of_trees = 5)
      predict(f, newdata = d$X_test)
    },
    predict_shrinkagetrees_interval_censored = function(baseline = FALSE) {
      d <- refactor_data(162)
      f <- fit(ShrinkageTrees, left_time = d$left, right_time = d$right,
               X_train = d$X, outcome_type = "interval-censored",
               prior_type = "horseshoe", number_of_trees = 5)
      predict(f, newdata = d$X_test)
    },
    predict_causalhorseforest_right_censored = function(baseline = FALSE) {
      d <- refactor_data(163)
      f <- fit(CausalHorseForest, y = d$time, status = d$status,
               X_train_control = d$X, X_train_treat = d$X,
               treatment_indicator_train = d$trt,
               outcome_type = "right-censored", number_of_trees = 5)
      predict(f, newdata_control = d$X_test, newdata_treat = d$X_test,
              bayesian_bootstrap = FALSE)
    }
  )
}

# Reduce a fitted or prediction object to the numbers worth comparing: every
# numeric top-level element, plus the preprocessing constants and prior scales
# stored for predict().
#
# Left out on purpose, because the refactor changed them by design:
#   * forestwide_shrinkage*   output of the removed "horseshoe_fw" prior
#   * zero-length elements    the old continuous causal fits carried two empty
#                             `*_predictions_sample` vectors
#   * args$a_dirichlet etc.   HorseTrees() used to record placeholder values
refactor_extract <- function(obj) {
  keep <- vapply(obj, function(x) is.numeric(x) && length(x) > 0, logical(1))
  out  <- unclass(obj)[keep]
  out  <- out[!grepl("^forestwide_shrinkage", names(out))]

  pre <- obj$preprocess
  if (!is.null(pre)) {
    out$meta_y_mean    <- pre$y_mean
    out$meta_sigma_hat <- pre$sigma_hat
  }

  a <- obj$args
  if (!is.null(a)) {
    out$meta_lambda <- a$lambda
    out$meta_scales <- c(a$param1, a$param2,
                         a$control$param1, a$control$param2,
                         a$treat$param1, a$treat$param2)
  }

  out[order(names(out))]
}

# Run one scenario quietly. Warnings are expected in places (for example a
# binary-looking continuous outcome) and are not what is being compared.
#
# The random number generator is pinned to R's defaults for the duration of
# the run. set.seed() alone is not enough: it fixes the seed, not the
# generator, so a session that has switched generator (for example to
# "L'Ecuyer-CMRG" for parallel work) would draw different data from the same
# seed and every scenario would differ from the baseline.
refactor_run <- function(scenario, baseline = FALSE) {
  old_kind <- RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  on.exit(suppressWarnings(RNGkind(old_kind[1], old_kind[2], old_kind[3])),
          add = TRUE)

  refactor_extract(suppressWarnings(scenario(baseline = baseline)))
}
