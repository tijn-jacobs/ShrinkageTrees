NewCausalShrinkageForest <- function(
  fit,
  call,
  outcome_type,
  timescale,
  n_train,
  p_control,
  p_treat,
  n_test,
  number_of_trees_control,
  number_of_trees_treat,
  N_post,
  N_burn,
  store_posterior_sample,
  sigma_hat,
  sigma_known,
  y_mean,
  prior_type_control_user,
  prior_type_control_cpp,
  prior_type_treat_user,
  prior_type_treat_cpp,
  dirichlet_bool_control = FALSE,
  dirichlet_bool_treat = FALSE
) {

  meta <- list(
    call = call,
    outcome_type = outcome_type,
    timescale = timescale,

    data_info = list(
      n_train = n_train,
      n_test = n_test,
      p_control = p_control,
      p_treat = p_treat
    ),

    mcmc = list(
      number_of_trees_control = number_of_trees_control,
      number_of_trees_treat = number_of_trees_treat,
      N_post = N_post,
      N_burn = N_burn,
      store_posterior_sample = store_posterior_sample
    ),

    preprocess = list(
      sigma_hat = sigma_hat,
      sigma_known = sigma_known,
      y_mean = y_mean
    ),

    prior = list(
      control = list(
        prior_type_user = prior_type_control_user,
        prior_type_cpp = prior_type_control_cpp,
        dirichlet = dirichlet_bool_control
      ),
      treat = list(
        prior_type_user = prior_type_treat_user,
        prior_type_cpp = prior_type_treat_cpp,
        dirichlet = dirichlet_bool_treat
      )
    )
  )

  obj <- c(fit, meta)
  class(obj) <- "CausalShrinkageForest"
  obj
}