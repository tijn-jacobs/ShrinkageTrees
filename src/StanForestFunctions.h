#ifndef GUARD_StanForestFunctions_h
#define GUARD_StanForestFunctions_h

#include "StanTree.h"
#include "StanTreeFunctions.h"
#include "Info.h"
#include <algorithm>

// Build a CutpointMatrix for p predictors from n observations stored in x
// (column-major).  nc[v] controls the number of cutpoints for predictor v:
// if nc != nullptr, uniform grids are used; if nc == nullptr, unique observed
// values are used as cutpoints.
void MakeCutpoints(size_t p, size_t n, double* x, CutpointMatrix& cutpoints,
                   int* nc);

// Compute the probability of proposing a birth step given the current tree.
// splittable_leaves is populated with all leaf nodes that can be split.
double GetBirthProbability(StanTree& tree, CutpointMatrix& cutpoints,
                           PriorInfo& prior_info,
                           std::vector<StanTree*>& splittable_leaves);

// Compute observation counts, weight sums and weighted residual sums for the
// left and right partitions induced by splitting target_leaf on
// (split_var, cut_val).  Weights default to 1 when data_info.weights is null.
void GetSufficientStatistics(StanTree& tree, StanTree* target_leaf,
                             size_t split_var, size_t cut_val,
                             CutpointMatrix& cutpoints, DataInfo& data_info,
                             size_t& left_count,
                             double& left_weight_sum,
                             double& left_weighted_residual_sum,
                             size_t& right_count,
                             double& right_weight_sum,
                             double& right_weighted_residual_sum);

// Compute observation counts, weight sums and weighted residual sums for an
// existing left/right leaf pair (used during a death proposal).
void GetSufficientStatistics(StanTree& tree, StanTree* left_leaf,
                             StanTree* right_leaf,
                             CutpointMatrix& cutpoints, DataInfo& data_info,
                             size_t& left_count,
                             double& left_weight_sum,
                             double& left_weighted_residual_sum,
                             size_t& right_count,
                             double& right_weight_sum,
                             double& right_weighted_residual_sum);

// Compute sufficient statistics for every leaf in the tree in a single pass
// over the data.  Populates leaves, weight_sums, and weighted_residual_sums.
void GetAllLeafStatistics(StanTree& tree, CutpointMatrix& cutpoints,
                          DataInfo& data_info,
                          std::vector<StanTree*>& leaves,
                          std::vector<double>& weight_sums,
                          std::vector<double>& weighted_residual_sums);

// Weighted leaf-marginal log-likelihood given Sum w_i and Sum w_i r_i under a
// Gaussian model with noise std dev sigma and leaf prior std dev eta.
double LogLikelihood(double weight_sum, double weighted_residual_sum,
                     double sigma, double eta);

// Probability that a node at the given depth grows; 0 if no valid split
// variable is available.
double ProbabilityNodeGrows(StanTree* node, CutpointMatrix& cutpoints,
                            PriorInfo& prior_info);

// Draw new step heights for all leaves of the tree from their posterior.
void DrawAllLeafMeans(StanTree& tree, CutpointMatrix& cutpoints,
                      DataInfo& data_info, PriorInfo& prior_info,
                      double sigma, Random& random);

// Generate a birth proposal: draw the target leaf, split variable, cut
// value, and compute the Metropolis ratio component (log_proposal_ratio).
void BirthProposal(StanTree& tree, CutpointMatrix& cutpoints,
                   PriorInfo& prior_info,
                   std::vector<StanTree*>& splittable_leaves,
                   double& prob_birth, StanTree*& target_leaf,
                   size_t& split_var, size_t& cut_val,
                   double& log_proposal_ratio,
                   std::vector<size_t>& variable_split_counts,
                   std::vector<double>& split_probabilities,
                   bool use_augmentation, Random& random);

// Generate a death proposal: draw the nog node to collapse and compute
// the Metropolis ratio component (log_proposal_ratio).
void DeathProposal(StanTree& tree, CutpointMatrix& cutpoints,
                   PriorInfo& prior_info,
                   std::vector<StanTree*>& splittable_leaves,
                   double& prob_birth, StanTree*& nog_node,
                   double& log_proposal_ratio, Random& random);

// Draw a single leaf mean from its Gaussian posterior.
double DrawLeafMean(double weight_sum, double weighted_residual_sum,
                    double eta, double sigma, Random& random);

// Draw the vector of splitting probabilities from a Dirichlet posterior
// (Linero, 2018).  Updates log_split_probabilities in-place.
void DrawSplitProbabilities(std::vector<size_t>& variable_split_counts,
                            std::vector<double>& log_split_probabilities,
                            double& dart_theta, Random& random);

// Draw the Dirichlet sparsity parameter theta from a grid approximation
// to its posterior (Linero, 2018).
void DrawSparsityParameter(bool fixed_theta, double& dart_theta,
                           std::vector<double>& log_split_probabilities,
                           double dart_a, double dart_b, double dart_rho,
                           Random& random);

#endif

