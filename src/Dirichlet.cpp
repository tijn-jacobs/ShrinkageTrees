// MIT License
// (c) Tijn Jacobs, 2025

#include "Dirichlet.h"


//
// Compute log-sum-exp in a numerically stable way.
// Used for normalizing log weights.
//
double LogSumExp(const vector<double>& values) {
  // Find maximum element for numerical stability
  double max_value = values[0];
  for (double x : values)
    if (x > max_value) max_value = x;

  // Compute sum of exponentials after centering by max_value
  double sum_exp = 0.0;
  for (double x : values)
    sum_exp += exp(x - max_value);

  // Return stabilized log-sum-exp
  return max_value + log(sum_exp);
}


//
// Draw log(s_j) from a Dirichlet(dirichlet_shape) distribution.
// Each component is sampled as gamma(shape_j, 1),
// then normalized to sum to 1, and log-transformed.
//
vector<double> LogDirichlet(const vector<double>& dirichlet_shape,
                            Random& random) {
  size_t p = dirichlet_shape.size();
  vector<double> samples(p);
  double sum_samples = 0.0;

  // Draw independent gamma variables
  for (size_t j = 0; j < p; ++j) {
    samples[j] = random.gamma(dirichlet_shape[j], 1.0);
    sum_samples += samples[j];
  }

  // Normalize and take logs
  for (size_t j = 0; j < p; ++j)
    samples[j] = log(samples[j] / sum_samples);

  return samples;
}


//
// DrawSplitProbs:
// Step 1 of the DART update.
// Given the number of times each covariate was used for splitting
// (variable_inclusion_count), sample a new vector of log splitting probabilities
// from the posterior Dirichlet distribution.
//
//   s | alpha, counts ~ Dirichlet(alpha/p + n_1, ..., alpha/p + n_p)
//
void DrawSplitProbs(vector<size_t>& variable_inclusion_count,
                    vector<double>& log_split_probs,
                    double& alpha_dirichlet,
                    Random& random) {
  
  // Compute number of variables
  size_t p = variable_inclusion_count.size();
  vector<double> alpha_post(p);

  // Compute posterior Dirichlet parameters
  for (size_t j = 0; j < p; ++j)
    alpha_post[j] = alpha_dirichlet / static_cast<double>(p)
                    + static_cast<double>(variable_inclusion_count[j]);

  // Sample log probabilities using Dirichlet draw
  log_split_probs = LogDirichlet(alpha_post, random);
}


//
// DrawDirichletAlpha:
// Step 2 of the DART update.
//
// Updates the concentration parameter alpha_dirichlet using
// a grid-based Gibbs sampler. The method follows the Bayesian
// sparsity scheme described by Linero (2018), but rewritten
// for ShrinkageTrees conventions.
//
// Model setup:
//   s ~ Dirichlet(alpha/p, ..., alpha/p)
//   lambda = alpha / (alpha + rho) ~ Beta(a, b)
//
// The posterior of alpha_dirichlet is approximated over a grid
// of lambda values in (0, 1).
//
void DrawDirichletAlpha(bool const_alpha,
                        double& alpha_dirichlet,
                        const vector<double>& log_split_probs,
                        vector<double>& variable_inclusion_prob,
                        double a_dirichlet,
                        double b_dirichlet,
                        double rho_dirichlet,
                        Random& random) {

  // Skip update when alpha is held constant
  if (const_alpha) return;

  const size_t num_vars = log_split_probs.size(); // ✅ define once here

  // Compute sufficient statistic: sum_j log(s_j)
  double log_sum_probs = 0.0;
  for (size_t j = 0; j < num_vars; ++j)
    log_sum_probs += log_split_probs[j];

  // Define lambda grid for the Beta(a, b) prior
  const size_t grid_points = 1000;
  vector<double> lambda_values(grid_points);
  vector<double> alpha_values(grid_points);
  vector<double> log_posterior(grid_points);

  // Evaluate the unnormalized log posterior on the grid
  for (size_t k = 0; k < grid_points; ++k) {
    // Map lambda in (0, 1) to alpha = (lambda * rho) / (1 - lambda)
    lambda_values[k] = static_cast<double>(k + 1) / (grid_points + 1.0);
    alpha_values[k] = (lambda_values[k] * rho_dirichlet) /
                      (1.0 - lambda_values[k]);

    // Log-likelihood term from the Dirichlet model
    double log_likelihood = lgamma(alpha_values[k])
      - num_vars * lgamma(alpha_values[k] / static_cast<double>(num_vars))
      + (alpha_values[k] / static_cast<double>(num_vars)) * log_sum_probs;

    // Log-prior term from Beta(a, b) on lambda
    double log_prior =
      (a_dirichlet - 1.0) * log(lambda_values[k]) +
      (b_dirichlet - 1.0) * log(1.0 - lambda_values[k]);

    log_posterior[k] = log_likelihood + log_prior;
  }

  // Normalize log weights using the stable log-sum-exp trick
  double log_norm = LogSumExp(log_posterior);
  for (size_t k = 0; k < grid_points; ++k)
    log_posterior[k] = exp(log_posterior[k] - log_norm);

  // Sample new alpha_dirichlet from the discrete grid posterior
  random.SetInclusionWeights(log_posterior);
  alpha_dirichlet = alpha_values[random.discrete()];

  // Back-transform log probabilities to normal scale
  variable_inclusion_prob.resize(num_vars);
  for (size_t j = 0; j < num_vars; ++j)
    variable_inclusion_prob[j] = exp(log_split_probs[j]);
}
