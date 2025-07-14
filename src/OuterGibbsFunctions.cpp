#include "OuterGibbsFunctions.h"

void UpdateForestwideShrinkage(
  std::vector<Tree>* trees,
  double& forestwide_shrinkage,
  double& forestwide_auxiliary,
  double alpha_fw,
  Rcpp::NumericVector& store_shrinkage,
  size_t i,
  size_t N_burn,
  Random& random
) {

  // Update auxiliary variable
  forestwide_auxiliary = random.inv_gamma(
    1.0, 1.0 / (alpha_fw * alpha_fw) + 1.0 / forestwide_shrinkage
  );

  // Compute sum of shrinkage contributions over all leaf nodes
  double sum_fw_shrink = 0.0;
  size_t n_leafs = 0;

  for (Tree& tree : *trees) {
    std::vector<Tree*> leaf_vector;
    tree.CollectLeafs(leaf_vector);

    for (Tree* leaf : leaf_vector) {
      double h = leaf->GetParameters(0);
      double lambda = leaf->GetParameters(1);
      sum_fw_shrink += (h * h) / lambda;
      n_leafs++;
    }
  }

  // Draw new forest-wide shrinkage parameter
  forestwide_shrinkage = random.inv_gamma(
    0.5 * (n_leafs + 1),
    0.5 * sum_fw_shrink + 1.0 / forestwide_auxiliary
  );

  // Store new shrinkage in the global parameter slot of each tree
  for (Tree& tree : *trees) {
    tree.SetGlobalParameters(0, forestwide_shrinkage);
    // Optionally also: tree.SetGlobalParameters(1, forestwide_auxiliary);
  }

  // Store result in posterior vector if we're past burn-in
  if (i >= N_burn) {
    store_shrinkage(i - N_burn) = forestwide_shrinkage;
  }
}

void AugmentCensoredObservations(
  bool is_survival,
  double* event_time, 
  const double* observed_time,
  const double* status_indicator,
  const double* predicted_time,
  const double& sigma,
  const size_t& n, 
  Random& random
) {

  if (!is_survival) return;

  // Declare variables
  double prediction, x_i, temporary, U;                     
  
  // Loop over all observations 
  for (size_t i = 0; i < n; i++) {

    // Check if the event time is censored
    if (status_indicator[i] == 0) {

      // Sample from the truncated normal distribution
      U = random.uniform();
      prediction = predicted_time[i];
      x_i = (observed_time[i] - prediction)/sigma;  

      // x_i is often way too large on the first few iterations !!!
      // if the normalized difference is greater than 4, set x_i = 4.0
      if(x_i > 4)  x_i = 4.0;

      temporary = U + (1.0 - U) * standard_normal_cdf(x_i);

      event_time[i] = prediction + sigma * standard_normal_quantile(temporary);
    }
  }
}

void UpdateSigma(
  bool sigma_known,
  double& sigma,
  Rcpp::NumericVector& store_sigma,
  const size_t i,
  const double* y,
  const size_t n,
  const double* prediction,
  const double nu,
  const double lambda,
  Random& random
) {

  if (sigma_known) return;

  double residual_sum2 = 0.0;
  double residual_temp;

  for (size_t k = 0; k < n; k++) {
    residual_temp = y[k] - prediction[k];
    residual_sum2 += residual_temp * residual_temp;
  }
  
  sigma = std::sqrt((nu * lambda + residual_sum2) / random.chi_square(n + nu));

  if (sigma >= 10) sigma = 10;

  store_sigma[i] = sigma;

  if (std::isnan(sigma) || std::isinf(sigma)) {
    Rcpp::stop("Sigma became invalid (NaN or Inf) during MCMC. Stopping execution.");
  }
}

void UpdateForestwideShrinkage(
  string prior_type,
  std::vector<Tree>* all_trees,
  Random& random,
  double& forestwide_auxiliary,
  double& forestwide_shrinkage,
  double alpha_fw
){

  if (prior_type != "horseshoe_fw") return;

  double sum_fw_shrink = 0.0;
  size_t n_leafs = 0;

  forestwide_auxiliary = random.inv_gamma(
    1.0,
    1.0 / (alpha_fw * alpha_fw) + 1.0 / forestwide_shrinkage
  );

  for (Tree& tree : *all_trees) {
    std::vector<Tree*> leaf_vector;
    tree.CollectLeafs(leaf_vector);

    for (Tree* leaf : leaf_vector) {
      double h = leaf->GetParameters(0);
      double lambda = leaf->GetParameters(1);
      sum_fw_shrink += (h * h) / lambda;
      n_leafs++;
    }
  }

  forestwide_shrinkage = random.inv_gamma(
    0.5 * (n_leafs + 1),
    0.5 * sum_fw_shrink + 1.0 / forestwide_auxiliary
  );

  for (Tree& tree : *all_trees) {
    tree.SetGlobalParameters(0, forestwide_shrinkage);
  }
}


