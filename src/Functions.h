// Summary: This header file defines functions for manipulating decision trees in a 
// Bayesian context. It includes functionality for identifying splittable nodes, 
// computing sufficient statistics, calculating GROW/PRUNE probabilities, and 
// updating tree parameters during MCMC sampling. Functions for printing tree 
// structures with node-specific information are also included.


#ifndef GUARD_Functions_H
#define GUARD_Functions_H

#include "Tree.h"

// Checks if a bottom node has variables that can be split on.
bool Splittable(Tree& leaf_node, Cutpoints& cutpoints);

// Collects all leaf nodes that can be split on.
void CollectSplittableLeaves(Tree& tree, std::vector<Tree*>& leaf_vector, 
                            Cutpoints& cutpoints);

// Finds variables that a node can split on and stores their indices in split_var.
void GetSplittableVariables(Tree& leaf, Cutpoints& cutpoints, 
                            std::vector<size_t>& split_var);

// Computes the probability of a GROW move with splittable_nodes containing 
// all valid bottom nodes.
double GrowProbability(Tree& tree, Cutpoints& cutpoints, 
                       TreePrior& tree_prior, 
                       std::vector<Tree*>& splittable_nodes);

// Computes the count (n) and sum of residuals (∑ y_i) for left and right 
// partitions given a split variable and cutpoint.
void SufficientStatistics(Tree& tree, Tree* target_node, size_t split_var, 
                          size_t cut_val, Cutpoints& cutpoints, Data& data, 
                          size_t& left_count, double& left_sum, 
                          size_t& right_count, double& right_sum);

// Computes the log-posterior likelihood of the data given the tree structure 
// and updated variance. Only used in non-reversible algorithm.
double LogPostLikelihood(size_t observation_count, double sum_residuals, 
                         double residual_std_dev, double prior_variance);

// Computes the probability that a node will grow; returns 0 if no valid 
// variables are available.
double ProbNodeGrows(Tree& node, Cutpoints& cutpoints, TreePrior& tree_prior);

// Computes the count (n) and sum of residuals (∑ y_i) for the left and right 
// child nodes of a split.
void SufficientStatistics(Tree& tree, Tree* left_child, Tree* right_child, 
                          Cutpoints& cutpoints, Data& data, size_t& left_count, 
                          double& left_sum, size_t& right_count, 
                          double& right_sum);

// Computes sufficient statistics for all bottom nodes in the tree.
void SufficientStatisticsAllLeaves(Tree& tree, Cutpoints& cutpoints, Data& data, 
                                   std::vector<Tree*>& bottom_nodes, 
                                   std::vector<size_t>& observation_count_vector, 
                                   std::vector<double>& residual_sum_vector);

// Draws mu (parameter) values for all bottom nodes in the tree.
void DrawMuAllLeaves(Tree& tree, Cutpoints& cutpoints, Data& data, 
                     TreePrior& tree_prior, double sigma, const double& omega,
                     ScaleMixture& scale_mixture, Random& random);

// Performs a full update of both global and local parameters for all bottom 
// nodes in the tree.
void FullUpdate(Tree& tree, 
                Cutpoints& cutpoints, 
                Data& data, 
                const double& sigma, 
                const double& omega, 
                ScaleMixture& scale_mixture, 
                Random& random);

// Draws a single mu value from the posterior distribution.
double DrawMuOneLeave(size_t n, 
                      double sum_residuals, 
                      double prior_variance, 
                      double sigma, 
                      Random& random);

// Prints the tree and includes the size of each leaf (bottom node) and the 
// residual sum.
void PrintTreeWithSizes(Tree& tree, Cutpoints& cutpoints, Data& data);

// Recursive function to print the tree structure with node sizes and residual sums.
void PrintTreeWithSizesRecursive(Tree* node, Cutpoints& cutpoints, Data& data, 
                                 std::map<const Tree*, size_t>& bottom_node_map, 
                                 std::vector<size_t>& observation_count_vector, 
                                 std::vector<double>& residual_sum_vector, 
                                 int depth);

// Computes the logarithm of the tree ratio for the PRUNE move.
double LogTreeRatio_PRUNE(Tree* node, const TreePrior& tree_prior, 
                          Cutpoints& cutpoints);

// Computes the logarithm of the tree ratio for the GROW move.
double LogTreeRatio_GROW(size_t depth, const TreePrior& tree_prior);

// Computes the logarithm of the move ratio for the GROW move.
double LogMoveRatio(size_t number_of_nogs, size_t number_of_leaves, 
                    double prune_prob, double grow_prob);

#endif // GUARD_Functions_H