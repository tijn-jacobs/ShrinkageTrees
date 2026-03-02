NewShrinkageTrees <- function(
  fit,
  call,
  outcome_type,
  timescale,
  prior_type_user,
  prior_type_cpp,
  n_train,
  p_features,
  n_test,
  number_of_trees,
  N_post,
  N_burn,
  store_posterior_sample,
  sigma_hat,
  sigma_known,
  y_mean
) {

  meta <- list(
    call = call,
    outcome_type = outcome_type,
    timescale = timescale,
    prior = list(
      prior_type_user = prior_type_user,
      prior_type_cpp = prior_type_cpp
    ),
    data_info = list(
      n_train = n_train,
      p_features = p_features,
      n_test = n_test
    ),
    mcmc = list(
      number_of_trees = number_of_trees,
      N_post = N_post,
      N_burn = N_burn,
      store_posterior_sample = store_posterior_sample
    ),
    preprocess = list(
      sigma_hat = sigma_hat,
      sigma_known = sigma_known,
      y_mean = y_mean
    )
  )

  obj <- c(fit, meta)

  class(obj) <- "ShrinkageTrees"

  obj
}