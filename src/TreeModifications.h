/* MIT License

 * Copyright (c) 2024 Tijn Jacobs

 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:

 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.

 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

// This header file declares functions for non-reversible and reversible tree 
// modifications in a Bayesian Additive Regression Tree (BART) framework. 
// The functions handle tree growth, pruning, changes, and reversible jumps,
// managing the selection of split variables and cutpoints, as well as 
// computing necessary statistics and acceptance ratios for Metropolis-Hastings 
// sampling. 

#ifndef GUARD_TreeModifications_H
#define GUARD_TreeModifications_H

#include "Functions.h"

// Function to grow a tree at a good bottom node
bool Grow(Tree& tree, Cutpoints& cutpoints, Data& data, TreePrior& tree_prior, 
          const double& sigma, std::vector<size_t>& variable_inclusion_count, 
          std::vector<double>& variable_inclusion_prob, 
          ScaleMixture& scale_mixture, Random& random);

// Function to prune a selected NOG node from the tree
bool Prune(Tree& tree, Cutpoints& cutpoints, Data& data, TreePrior& tree_prior, 
           const double& sigma, std::vector<size_t>& variable_inclusion_count, 
           std::vector<double>& variable_inclusion_prob, 
           ScaleMixture& scale_mixture, Random& random);

// Function to change the split variable and cut value of a NOG node
bool Change(Tree& tree, Cutpoints& cutpoints, Data& data, TreePrior& tree_prior, 
            const double& sigma, std::vector<size_t>& variable_inclusion_count, 
            std::vector<double>& variable_inclusion_prob, 
            ScaleMixture& scale_mixture, Random& random);

// Function to perform a reversible jump grow operation
bool RJ_Grow(Tree& tree, Cutpoints& cutpoints, Data& data, TreePrior& tree_prior, 
             const double& sigma, std::vector<size_t>& variable_inclusion_count, 
             std::vector<double>& variable_inclusion_prob, 
             ScaleMixture& scale_mixture, Random& random);

// Function to perform a reversible jump prune operation
bool RJ_Prune(Tree& tree, Cutpoints& cutpoints, Data& data, TreePrior& tree_prior, 
              const double& sigma, std::vector<size_t>& variable_inclusion_count, 
              std::vector<double>& variable_inclusion_prob, 
              ScaleMixture& scale_mixture, Random& random);

// Function to perform a reversible jump change operation
bool RJ_Change(Tree& tree, Cutpoints& cutpoints, Data& data, TreePrior& tree_prior, 
               const double& sigma, std::vector<size_t>& variable_inclusion_count, 
               std::vector<double>& variable_inclusion_prob, 
               ScaleMixture& scale_mixture, Random& random);

// Function to perform reversible jumps with specific probabilities for each action
bool ReversibleJump(Tree& tree, Cutpoints& cutpoints, Data& data, 
                    TreePrior& tree_prior, const double& sigma, const double& omega, 
                    std::vector<size_t>& nodes, 
                    std::vector<double>& variable_inclusion_prob, 
                    const size_t& delayed_proposal,
                    bool reversible_jump, ScaleMixture& scale_mixture, 
                    Random& random);

#endif // GUARD_TreeModifications_H