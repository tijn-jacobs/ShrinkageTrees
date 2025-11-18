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
  for (size_t j = 0; j < p; ++j)  samples[j] = log(samples[j] / sum_samples);

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
                    vector<double>& log_vip,
                    double& alpha_dirichlet,
                    Random& random) {
  

// cout << "Drawing split probs with alpha_dirichlet = " << alpha_dirichlet << "\n";
  
// Compute number of variables
  size_t p = variable_inclusion_count.size();
  vector<double> alpha_post(p);

  // Compute posterior Dirichlet parameters
  for (size_t j = 0; j < p; ++j)
    alpha_post[j] = alpha_dirichlet / static_cast<double>(p)
                    + static_cast<double>(variable_inclusion_count[j]);

  // Sample log probabilities using Dirichlet draw
  log_vip = LogDirichlet(alpha_post, random);
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

/*
void DrawDirichletAlpha(bool const_alpha,
                        double& alpha_dirichlet,
                        vector<double>& log_vip,
                        double a_dirichlet,
                        double b_dirichlet,
                        double rho_dirichlet,
                        Random& random) {

  // Skip update when alpha is held constant
  if (const_alpha) return;

  cout << "Drawing Dirichlet alpha with current value = " << alpha_dirichlet << "\n";

  const size_t num_vars = log_vip.size(); // ✅ define once here

  // Compute sufficient statistic: sum_j log(s_j)
  double log_sum_probs = 0.0;
  for (size_t j = 0; j < num_vars; ++j)
    log_sum_probs += log_vip[j];

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
}
*/

double log_sum_exp(std::vector<double>& v){
    double mx=v[0],sm=0.;
    for(size_t i=0;i<v.size();i++) if(v[i]>mx) mx=v[i];
    for(size_t i=0;i<v.size();i++){
      sm += exp(v[i]-mx);
    }
    return mx+log(sm);
}

void DrawDirichletAlpha(bool const_theta, double& theta, std::vector<double>& lpv,
		 double a, double b, double rho, Random& random){

  // Draw sparsity parameter theta_0 (Linero calls it alpha); see Linero, 2018
  // theta / (theta + rho ) ~ Beta(a,b)
  // Set (a=0.5, b=1) for sparsity
  // Set (a=1, b=1) for non-sparsity
  // rho = p usually, but making rho < p increases sparsity
  if(const_theta) return;

  size_t p=lpv.size();
  double sumlpv=0.,lse;
    
  std::vector<double> lambda_g (1000,0.);
  std::vector<double> theta_g (1000,0.);
  std::vector<double> lwt_g (1000,0.);
  for(size_t j=0;j<p;j++) sumlpv+=lpv[j];

  for(size_t k=0;k<1000;k++){
    lambda_g[k]=(double)(k+1)/1001.;
    theta_g[k]=(lambda_g[k]*rho)/(1.-lambda_g[k]);
    double theta_log_lik=lgamma(theta_g[k])-(double)p*lgamma(theta_g[k]/(double)p)+(theta_g[k]/(double)p)*sumlpv;
    double beta_log_prior=(a-1.)*log(lambda_g[k])+(b-1.)*log(1.-lambda_g[k]);
//      cout << "SLP: " << sumlogpv << "\nTLL: " << theta_log_lik << "\nBLP: " << beta_log_prior << '\n';
    lwt_g[k]=theta_log_lik+beta_log_prior;      
  }
  
  lse=log_sum_exp(lwt_g);

  for(size_t k=0;k<1000;k++) {
    lwt_g[k]=exp(lwt_g[k]-lse);
//      cout << "LWT: " << lwt_g[k] << '\n';
  }

  random.SetInclusionWeights(lwt_g);    
  theta=theta_g[random.discrete()];
}