#include "HorseTrees.h"
#include "Timing.h"


// [[Rcpp::export]]
Rcpp::List HorseTrees_cpp( 
  SEXP nSEXP,     
  SEXP pSEXP,     
  SEXP n_testSEXP,    
  SEXP X_trainSEXP,    
  SEXP ySEXP,  
  SEXP status_indicatorSEXP, 
  SEXP is_survivalSEXP,     
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
  SEXP sigmaSEXP,
  SEXP sigma_knownSEXP,     
  SEXP omegaSEXP,
  SEXP param1SEXP,
  SEXP param2SEXP,
  SEXP prior_typeSEXP,
  SEXP reversibleSEXP,
  SEXP store_parametersSEXP,
  SEXP store_posterior_sampleSEXP,
  SEXP n1SEXP,     
  SEXP n2SEXP,
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
  double sigma = Rcpp::as<double>(sigmaSEXP);
  bool sigma_known = Rcpp::as<bool>(sigma_knownSEXP);
  double omega = Rcpp::as<double>(omegaSEXP);
  bool reversible = Rcpp::as<bool>(reversibleSEXP);
  
  double param1 = Rcpp::as<double>(param1SEXP);   
  double param2 = Rcpp::as<double>(param2SEXP);
  string prior_type = Rcpp::as<string>(prior_typeSEXP);
  bool store_parameters =  Rcpp::as<bool>(store_parametersSEXP);

  unsigned int n1 = Rcpp::as<unsigned int>(n1SEXP);   
  unsigned int n2 = Rcpp::as<unsigned int>(n2SEXP);  
  bool print_progress = Rcpp::as<bool>(verboseSEXP);

  // Objects for is_survival data
  bool is_survival = Rcpp::as<bool>(is_survivalSEXP);
  Rcpp::NumericVector status_indicator_vector(status_indicatorSEXP);   
  double* status_indicator = &status_indicator_vector[0];   
  std::vector<double> y_observed_vector(y_vector.begin(), y_vector.end());  // Create an independent copy of y
  double* y_observed = y_observed_vector.data(); // Get pointer to use like an array

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
  
  std::vector<size_t> cumulative_inclusion_count(p, 0);
  std::vector<std::vector<double>> variable_inclusion_prob;
  
  // Initialize the C-based random number generator
  arn random(n1, n2);

  // Map the input to the corresponding PriorType
  PriorType prior;
  
  if (prior_type == "horseshoe") {
    prior = PriorType::Horseshoe;
  } else if (prior_type == "fixed" || prior_type == "standard") {
    prior = PriorType::FixedVariance;
  } else if (prior_type == "halfcauchy") {
    prior = PriorType::HalfCauchy;
  } else if (prior_type == "horseshoe_fw") {
    prior = PriorType::Horseshoe_fw;
  } else {
    Rcpp::stop("Invalid prior type provided. Choose one of: 'horseshoe', 'fixed', 'halfcauchy', 'horseshoe_fw', 'standard'.");
  }

  // What happens with prior if we use 'standard', i.e., non-reversible-jump BART?
  // It will not be used. So can we set it to "nothing"?
  
  // Initialize the scale mixture prior on the step heights in the leaves
  ScaleMixture scale_mixture(prior, param1, param2);
  
  // Build the forest
  Forest forest(number_of_trees);
  forest.SetTreePrior(base, power, param1, p_grow, p_prune); // In case of NON-RJ; param1 = step height variance
  forest.SetUpForest(p, n, X_train, y, nullptr, omega);
  
  for (size_t i = 0; i < n; i++) train_predictions_mean[i] = 0.0;
  for (size_t i = 0; i < n_test; i++) test_predictions_mean[i] = 0.0;
  
  //-----------------------------------------------------------
  
  // Temporary storage
  double* testpred = (n_test) ? new double[n_test] : nullptr;
  
  //-----------------------------------------------------------
  // MCMC
  if(print_progress) Rcpp::Rcout << "\nProgress of the MCMC sampler:\n\n";

  
  size_t total = N_post + N_burn;
  int barWidth = 70;
  
  bool* accepted = new bool[number_of_trees];
  double sum_accept = 0;
  double acceptance_ratio;

  int max_stored_leaves = 1;
  if(store_parameters)  max_stored_leaves = 20;
  std::vector<Tree>* all_trees = forest.GetTreesPointer();

  Rcpp::NumericMatrix store_global_parameters;
  Rcpp::NumericMatrix store_local_parameters;
  Rcpp::IntegerMatrix store_local_indices;

  if (store_parameters) {
    store_global_parameters = Rcpp::NumericMatrix(N_post, number_of_trees);
    store_local_parameters = Rcpp::NumericMatrix(N_post, max_stored_leaves * number_of_trees);
    store_local_indices = Rcpp::IntegerMatrix(N_post, max_stored_leaves * number_of_trees);
    
    // Initialize with -2 for debugging (optional)
    store_global_parameters.fill(-2);
    store_local_parameters.fill(-2);
    store_local_indices.fill(-2);
  }

  // Forest wide Horseshoe shrinkage parameter and aux.
  double forestwide_shrinkage = 1.0;
  double forestwide_auxiliary = 1.0;
  Rcpp::NumericVector store_forestwide_shrinkage;
  
  if (prior_type == "horseshoe_fw") {
    store_forestwide_shrinkage = Rcpp::NumericVector(N_post);
    store_forestwide_shrinkage.fill(0);
  }  

  time_t time_stamp;
  int time1 = time(&time_stamp);
  
  for (size_t i = 0; i < total; i++) {
    
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


    {
    ScopedTimer t("UpdateForest");
    // Update the forest (outer Gibbs step)
    forest.UpdateForest(sigma, scale_mixture, reversible, delayed_proposal, random, accepted);
    }
  
    // Update fores wide shrinkage parameter (outer Gibbs step0
    if (prior_type == "horseshoe_fw") {
      UpdateForestwideShrinkage(
        all_trees,
        forestwide_shrinkage,
        forestwide_auxiliary,
        param2,
        store_forestwide_shrinkage,
        i,
        N_burn,
        random
      );
    }
  
    {
    ScopedTimer t("UpdateSigma");
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
    }

    {
    ScopedTimer t("AugmentCensoredObservations");
    // Augment the censored data
    AugmentCensoredObservations(is_survival, y, y_observed, status_indicator, forest.GetPredictions(), sigma, n, random);
    }

    {
    ScopedTimer t("StoreParametersAuxilliary");
    // Save the (averages for now) of the leaf node parameters
    if (store_parameters && i >= N_burn) {
      size_t tree_counter = 0;

      for (Tree& tree : *all_trees) {
        // Collect all leaf nodes
        std::vector<Tree*> leaf_vector;
        tree.CollectLeaves(leaf_vector);

        // Iterate through the collected leaf nodes
        for (size_t leaf_index = 0; leaf_index < leaf_vector.size(); ++leaf_index) {
          Tree* leaf = leaf_vector[leaf_index];
              
          int split_var_int = - 1;

          // Get the single parameter for the current leaf
          double parameter = leaf->GetParameters(1);  // Assuming GetParameters() returns a single double
          if (leaf->GetParent()) {
            split_var_int = static_cast<int>(leaf->GetParent()->GetSplitVar());
          } 

          // Compute the column index for this leaf
          size_t col_index = tree_counter * max_stored_leaves + leaf_index;

          // Ensure col_index does not exceed matrix bounds
          if (col_index < static_cast<size_t>(store_local_parameters.ncol())) {
              store_local_parameters(i - N_burn, col_index) = parameter;
              store_local_indices(i - N_burn, col_index) = split_var_int;
          } else {
              Rcpp::stop("Index out of bounds while storing leaf parameters.");
          }
        }

        // Retrieve the global parameters from the root node and store it
        store_global_parameters(i - N_burn, tree_counter) = tree.GetGlobalParameters(0);

        // Keep track of the trees
        tree_counter++;
      }
    }
    }

    {
    ScopedTimer t("StoreParameters");
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

      // Accumulate variable inclusion counts
      for (size_t var = 0; var < p; ++var) {
        cumulative_inclusion_count[var] += forest.GetVariableInclusionCount()[var];
      }

      // Predict test set
      if (n_test > 0) {
        forest.Predict(p, n_test, X_test, testpred);
      }

      // Store posterior sample of test predictions
      if (store_posterior_sample && n_test > 0) {
        for (size_t k = 0; k < n_test; k++) {
          test_predictions_sample(i - N_burn, k) = testpred[k];
        }
      }

      // Store posterior mean of test predictions
      if (n_test > 0) {
        for (size_t k = 0; k < n_test; k++) {
          test_predictions_mean[k] += testpred[k];
        }
      }

      // Track acceptance
      for (size_t j = 0; j < number_of_trees; j++) {
        sum_accept += accepted[j];
      }
    }
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
  if (store_parameters) {
    results["global_shrinkage_parameters"] = store_global_parameters;
    results["local_shrinkage_parameters"] = store_local_parameters;
    results["local_splitting_variables"] = store_local_indices;
  }
  if (prior_type == "horseshoe_fw") {
    results["forestwide_shrinkage"] = store_forestwide_shrinkage;
  } 
  
  if (testpred) delete[] testpred;
  delete[] accepted;

  return results;
}