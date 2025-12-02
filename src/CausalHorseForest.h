#ifndef CAUSAL_HORSETREES_H
 #define CAUSAL_HORSETREES_H
 
 #include "HorseTrees.h"
 #include "probit-HorseTrees.h"
 #include "OuterGibbsFunctions.h"
 
 
  Rcpp::List CausalHorseForest_cpp(
    SEXP nSEXP, SEXP p_treatSEXP, SEXP p_controlSEXP, SEXP X_train_treatSEXP,
    SEXP X_train_controlSEXP, SEXP ySEXP, SEXP status_indicatorSEXP, SEXP is_survivalSEXP,
    SEXP treatment_indicatorSEXP,
    SEXP n_testSEXP, SEXP X_test_controlSEXP, SEXP X_test_treatSEXP,
    SEXP treatment_indicator_testSEXP, 
    SEXP no_trees_treatSEXP, SEXP power_treatSEXP, SEXP base_treatSEXP,
    SEXP p_grow_treatSEXP, SEXP p_prune_treatSEXP, SEXP omega_treatSEXP,
    SEXP prior_type_treatSEXP, SEXP param1_treatSEXP, SEXP param2_treatSEXP,
    SEXP reversible_treatSEXP, SEXP no_trees_controlSEXP,
    SEXP power_controlSEXP, SEXP base_controlSEXP, SEXP p_grow_controlSEXP,
    SEXP p_prune_controlSEXP, SEXP omega_controlSEXP, SEXP prior_type_controlSEXP,
    SEXP param1_controlSEXP, SEXP param2_controlSEXP, SEXP reversible_controlSEXP,
    SEXP sigma_knownSEXP, SEXP sigmaSEXP, SEXP lambdaSEXP,
    SEXP nuSEXP, SEXP N_postSEXP, SEXP N_burnSEXP, SEXP delayed_proposalSEXP,
    SEXP store_parametersSEXP, SEXP max_stored_leavesSEXP,
    SEXP store_posterior_sampleSEXP, SEXP verboseSEXP
  );
 
 #endif // CAUSAL_HORSETREES_H
 