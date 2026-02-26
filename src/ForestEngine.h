#ifndef GUARD_ForestEngine_h
#define GUARD_ForestEngine_h


#include "StanForest.h"


// This class handles the translation from the StanForest implementation for 
// the standard BART and dirichlet BART models, and the Forest implementation for
// the ShrinkageTrees models with reversible jump moves.

enum class ForestEngineType {
  forest_type,
  stan_forest_type
};

struct ForestEngine {

  ForestEngineType type;

  std::unique_ptr<Forest> forest;
  std::unique_ptr<StanForest> stan_forest;

  ForestEngine(ForestEngineType type, size_t m) : type(type) {
    if (type == ForestEngineType::forest_type) {
      forest = std::make_unique<Forest>(m);
    } else {
      stan_forest = std::make_unique<StanForest>(m);
    }
  }

  // Set-up the forest
  void SetTreePrior(double base, double power, double eta, double p_grow, 
                   double p_prune, double a_dirichlet, double b_dirichlet,
                   double rho_dirichlet, bool augment, bool dirichlet_bool,
                   double alpha_dirichlet) {
                      
    if (type == ForestEngineType::forest_type) {
      forest->SetTreePrior(base, power, eta, p_grow, p_prune);
    } else {
      stan_forest->setprior(base, power, eta);   
      stan_forest->setdart(a_dirichlet, b_dirichlet, rho_dirichlet, augment, dirichlet_bool, alpha_dirichlet);
    } // Usually, rho_dirichlet = p, and augment = true
  }

  void SetUpForest(size_t p, size_t n, double* X, double* augment_outcome, int* nc, 
                   double omega) { // the augment_outcome does not have to be augmented, BUT CAN BE!

    if (type == ForestEngineType::forest_type) {
      forest->SetUpForest(p, n, X, augment_outcome, nc, omega);
    } else {
      stan_forest->setdata(p, n, X, augment_outcome, nc);
    }
  }

  void StartDirichlet() {
    if (type == ForestEngineType::stan_forest_type) {
      stan_forest->startdart();
    }
  }

  // Update the forest
  void UpdateForest(double sigma,
              ScaleMixture& scale_mixture,
              bool reversible,
              size_t delayed,
              Random& rng,
              bool* accepted) 
  {
    if (type == ForestEngineType::forest_type) {
      forest->UpdateForest(
        sigma, scale_mixture, reversible, delayed, rng, accepted
      );
    } else {
      stan_forest->draw(sigma, rng, accepted);  // StanForest does not use scale mixture
      // optionally: fill accepted[] = ??? (Stan doesn't use it)
    }
  }

  // Predict new outcomes
  void Predict(size_t p, size_t n_test, double* X, double* out) {
    if (type == ForestEngineType::forest_type) {
      forest->Predict(p, n_test, X, out);
    } else {
      stan_forest->predict(p, n_test, X, out);
    }
  }

  // Accessors
  inline double GetPrediction(size_t i) {
    return (type == ForestEngineType::forest_type)
               ? forest->GetPrediction(i)
               : stan_forest->f(i);
  }

  inline double* GetPredictions() {
    return (type == ForestEngineType::forest_type)
               ? forest->GetPredictions()
               : stan_forest->GetAllFit();
  }

  std::vector<size_t>& GetVariableInclusionCount() {
    if (type == ForestEngineType::forest_type) {
      return forest->GetVariableInclusionCount();
    } else {
      return stan_forest->getnv();   // correct Stan version
    }
  }

  std::vector<double>& GetVariableInclusionProb() {
    if (type == ForestEngineType::forest_type) {
      return forest->GetVariableInclusionProb();
    } else {
      return stan_forest->getpv();   // correct Stan version
    }
  }

  void UpdateGlobalScaleParameters(string prior_type,
                                   double global_parameter,
                                   double& storage_eta, // Store the updated eta at this location
                                   Random& random) {
    if (type == ForestEngineType::forest_type) {
      return;
    } else {
      return stan_forest->UpdateGlobalScaleParameters(prior_type, global_parameter, storage_eta, random);   // correct Stan version
    }
  }
};

#endif
