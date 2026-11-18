#include "HorseTrees.h"
#include "ForestEngine.h"

// [[Rcpp::export]]
Rcpp::List HorseTrees_cpp( 
  SEXP nSEXP,     
  SEXP pSEXP,     
  SEXP n_testSEXP,    
  SEXP X_trainSEXP,    
  SEXP ySEXP,  
  SEXP status_indicatorSEXP,
  SEXP is_survivalSEXP,
  SEXP observed_left_timeSEXP,
  SEXP observed_right_timeSEXP,
  SEXP interval_censoring_indicatorSEXP,
  SEXP X_testSEXP,    
  SEXP number_of_treesSEXP,     
  SEXP N_postSEXP,     
  SEXP N_burnSEXP,  
  SEXP delayed_proposalSEXP, 
  SEXP powerSEXP,    
  SEXP baseSEXP,
  SEXP p_growSEXP,
  SEXP p_pruneSEXP,
  SEXP nuSEXP,     
  SEXP lambdaSEXP,
  SEXP dirichlet_boolSEXP,
  SEXP a_dirichletSEXP,
  SEXP b_dirichletSEXP, 
  SEXP rho_dirichletSEXP,     
  SEXP sigmaSEXP,
  SEXP sigma_knownSEXP,     
  SEXP omegaSEXP,
  SEXP param1SEXP,
  SEXP param2SEXP,
  SEXP prior_typeSEXP,
  SEXP reversibleSEXP,
  SEXP store_posterior_sampleSEXP,
  SEXP verboseSEXP
) {
  
  // Explicit conversion of SEXP to appropriate C++ types using Rcpp::as   
  size_t n = Rcpp::as<size_t>(nSEXP);   
  size_t p = Rcpp::as<size_t>(pSEXP);   
  size_t n_test = Rcpp::as<size_t>(n_testSEXP);    
  
  // Converting SEXP to Rcpp NumericVector for matrix-like data
  Rcpp::NumericVector X_train_vector(X_trainSEXP);   
  double* X_train = &X_train_vector[0];   
  
  Rcpp::NumericVector y_vector(ySEXP);   
  double* y = &y_vector[0];   
  
  Rcpp::NumericVector X_test_vector(X_testSEXP);   
  double* X_test = &X_test_vector[0];    
  
  size_t number_of_trees = Rcpp::as<size_t>(number_of_treesSEXP);    
  
  
  size_t N_post = Rcpp::as<size_t>(N_postSEXP);   
  size_t N_burn = Rcpp::as<size_t>(N_burnSEXP);   
  size_t delayed_proposal = Rcpp::as<size_t>(delayed_proposalSEXP);   
  double power = Rcpp::as<double>(powerSEXP);   
  double base = Rcpp::as<double>(baseSEXP);  
  double p_grow = Rcpp::as<double>(p_growSEXP);   
  double p_prune = Rcpp::as<double>(p_pruneSEXP);   
  double nu = Rcpp::as<double>(nuSEXP);   
  double lambda = Rcpp::as<double>(lambdaSEXP);  
  bool dirichlet_bool = Rcpp::as<bool>(dirichlet_boolSEXP);
  double a_dirichlet = Rcpp::as<double>(a_dirichletSEXP);
  double b_dirichlet = Rcpp::as<double>(b_dirichletSEXP); 
  double rho_dirichlet = Rcpp::as<double>(rho_dirichletSEXP);
  double sigma = Rcpp::as<double>(sigmaSEXP);
  bool sigma_known = Rcpp::as<bool>(sigma_knownSEXP);
  double omega = Rcpp::as<double>(omegaSEXP);
  bool reversible = Rcpp::as<bool>(reversibleSEXP);
  
  double param1 = Rcpp::as<double>(param1SEXP);   
  double param2 = Rcpp::as<double>(param2SEXP);
  string prior_type = Rcpp::as<string>(prior_typeSEXP);

  bool print_progress = Rcpp::as<bool>(verboseSEXP);

  // Objects for is_survival data
  bool is_survival = Rcpp::as<bool>(is_survivalSEXP);
  Rcpp::NumericVector status_indicator_vector(status_indicatorSEXP);
  double* status_indicator = &status_indicator_vector[0];

  // Objects for interval censoring
  Rcpp::NumericVector observed_left_time_vector(observed_left_timeSEXP);
  double* observed_left_time = &observed_left_time_vector[0];
  Rcpp::NumericVector observed_right_time_vector(observed_right_timeSEXP);
  double* observed_right_time = &observed_right_time_vector[0];
  Rcpp::NumericVector interval_censoring_indicator_vector(interval_censoring_indicatorSEXP);
  double* interval_censoring_indicator = &interval_censoring_indicator_vector[0];

  // Declare sigma storage based on whether sigma is known
  Rcpp::NumericVector store_sigma = sigma_known ? Rcpp::NumericVector::create(sigma) : Rcpp::NumericVector(N_post + N_burn);

  // Declare what has to be returned
  Rcpp::NumericVector train_predictions_mean(n);
  Rcpp::NumericVector test_predictions_mean(n_test);

  bool store_posterior_sample = Rcpp::as<bool>(store_posterior_sampleSEXP);

  Rcpp::NumericMatrix train_predictions_sample;
  Rcpp::NumericMatrix test_predictions_sample;

  if (store_posterior_sample) {
    train_predictions_sample = Rcpp::NumericMatrix(N_post, n);
    test_predictions_sample = Rcpp::NumericMatrix(N_post, n_test);
  }

  Rcpp::NumericVector store_alpha_dirichlet;
  Rcpp::NumericMatrix store_split_probs;
  Rcpp::NumericMatrix store_split_counts;

  double alpha_dirichlet = 1.0;
  if (dirichlet_bool) {
    store_alpha_dirichlet = Rcpp::NumericVector(N_post);
    store_split_probs     = Rcpp::NumericMatrix(N_post, p);
  }
  store_split_counts     = Rcpp::NumericMatrix(N_post, p);

  // Initialize the C-based random number generator
  RandomGenerator random;

  // Map the input to the corresponding PriorType
  PriorType prior;
  
  if (prior_type == "horseshoe") {
    prior = PriorType::Horseshoe;
  } else if (prior_type == "standard") {
    prior = PriorType::FixedVariance;
  } else if (prior_type == "halfcauchy") {
    prior = PriorType::HalfCauchy;
  } else if (prior_type == "dirichlet") {
    prior = PriorType::FixedVariance; // <- This should be an empty prior!!!
  } else {
    Rcpp::stop("Invalid prior type provided. Choose one of: 'horseshoe', 'halfcauchy', 'standard', 'dirichlet'.");
  }
  
  // Initialize the scale mixture prior on the step heights in the leaves
  ScaleMixture scale_mixture(prior, param1, param2);
  
  // Build the forest
  ForestEngineType engine_type;
  if (reversible) {
    engine_type = ForestEngineType::forest_type;
  } else {
    engine_type = ForestEngineType::stan_forest_type;
  }
  ForestEngine forest(engine_type, number_of_trees);
  forest.SetTreePrior(base, power, param1, p_grow, p_prune, a_dirichlet, b_dirichlet, rho_dirichlet, false, dirichlet_bool, alpha_dirichlet); // In case of NON-RJ; param1 = step height variance
  forest.SetUpForest(p, n, X_train, y, nullptr, omega);

  
  for (size_t i = 0; i < n; i++) train_predictions_mean[i] = 0.0;
  for (size_t i = 0; i < n_test; i++) test_predictions_mean[i] = 0.0;
  
  
  // Temporary storage
  double* test_predictions = (n_test) ? new double[n_test] : nullptr;
  
  // MCMC
  if(print_progress) Rcpp::Rcout << "\nProgress of the MCMC sampler:\n\n";

  
  size_t total = N_post + N_burn;
  int barWidth = 70;
  
  bool* accepted = new bool[number_of_trees];
  double sum_accept = 0;
  double acceptance_ratio;

  time_t time_stamp;
  int time1 = time(&time_stamp);
  
  for (size_t i = 0; i < total; i++) {

    if(i == (N_burn/2) && dirichlet_bool) {
      forest.StartDirichlet();
    }

    if(print_progress){
      // Progress bar
      float progress = static_cast<float>(i) / static_cast<float>(total);
      Rcpp::Rcout << "|";
      int pos = static_cast<int>(barWidth * progress);
      for (int j = 0; j < barWidth; ++j) {
        if (j < pos) Rcpp::Rcout << "=";
        else if (j == pos) Rcpp::Rcout << ">";
        else Rcpp::Rcout << " ";
      }
      Rcpp::Rcout << "| " << int(progress * 100.0) << " %\r";
      Rcpp::Rcout.flush();
    }

    // Update the forest (outer Gibbs step)
    forest.UpdateForest(
      sigma, 
      scale_mixture, 
      reversible, 
      delayed_proposal, 
      random, 
      accepted
    );
    
    // Update sigma (outer Gibbs step)
    UpdateSigma(
      sigma_known,
      sigma,
      store_sigma,
      i,
      y,
      n,
      forest.GetPredictions(),
      nu,
      lambda,
      random
    );


    // Augment the censored data (right-censored and interval-censored)
    AugmentCensoredObservations(
      is_survival,
      y,
      observed_left_time,
      status_indicator,
      observed_right_time,
      interval_censoring_indicator,
      forest.GetPredictions(),
      sigma,
      n,
      random
    );

    if (i >= N_burn) {
      // Store posterior mean of training predictions
      for (size_t k = 0; k < n; k++) {
        train_predictions_mean[k] += forest.GetPrediction(k);
      }

      // Store posterior sample of training predictions
      if (store_posterior_sample) {
        for (size_t k = 0; k < n; k++) {
          train_predictions_sample(i - N_burn, k) = forest.GetPrediction(k);
        }
      }

      // Predict test set
      if (n_test > 0) {
        forest.Predict(p, n_test, X_test, test_predictions);
      }

      // Store posterior sample of test predictions
      if (store_posterior_sample && n_test > 0) {
        for (size_t k = 0; k < n_test; k++) {
          test_predictions_sample(i - N_burn, k) = test_predictions[k];
        }
      }

      // Store posterior mean of test predictions
      if (n_test > 0) {
        for (size_t k = 0; k < n_test; k++) {
          test_predictions_mean[k] += test_predictions[k];
        }
      }

      // Track acceptance
      for (size_t j = 0; j < number_of_trees; j++) {
        sum_accept += accepted[j];
      }

      if (dirichlet_bool) {
        const std::vector<double>& probs = forest.GetVariableInclusionProb();
        for (size_t j = 0; j < p; ++j) {
          store_split_probs(i - N_burn, j) = probs[j];
        }
      }
      const std::vector<size_t>& counts = forest.GetVariableInclusionCount();
      for (size_t j = 0; j < p; ++j) store_split_counts(i - N_burn, j) = counts[j];
      
    }
  }

  for (size_t k = 0; k < n; k++) train_predictions_mean[k] /= N_post;
  for (size_t k = 0; k < n_test; k++) test_predictions_mean[k] /= N_post;

  // Calculate the acceptance ratio
  acceptance_ratio = sum_accept / (N_post * number_of_trees);

  int time2 = time(&time_stamp);
  
    
  if(print_progress) {
    // Final update of progress bar to ensure it ends at 100%
    Rcpp::Rcout << "|";
    for (int j = 0; j < barWidth; ++j) {
      Rcpp::Rcout << "=";
    }
    Rcpp::Rcout << "| 100 %\r";
    Rcpp::Rcout.flush();
    
    Rcpp::Rcout << "\n" << std::endl;  // Move to the next line after the progress bar is done.

    Rcpp::Rcout << "Mean acceptance ratio: " << acceptance_ratio << std::endl;
    Rcpp::Rcout << "\nDone in " << (time2 - time1) << " seconds.\n";
    Rcpp::Rcout << std::endl;
  }


  Rcpp::List results;
  results["sigma"] = store_sigma;
  results["test_predictions"] = test_predictions_mean;
  results["train_predictions"] = train_predictions_mean;
  results["acceptance_ratio"] = acceptance_ratio;
  if (store_posterior_sample) {
    results["test_predictions_sample"] = test_predictions_sample;
    results["train_predictions_sample"] = train_predictions_sample;
  } 
  if (dirichlet_bool) {
    results["alpha_dirichlet"] = store_alpha_dirichlet;
    results["split_probs"] = store_split_probs;
  }
  results["store_split_counts"] = store_split_counts;

  
  if (test_predictions) delete[] test_predictions;
  delete[] accepted;

  return results;
}
