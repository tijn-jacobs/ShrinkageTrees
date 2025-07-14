// Computational assumption: any node can be split. That is, there always exists a covariate that can be chosen to assign a split rule to.
// This boils down to having "enough" continuous covariates in the data. 

// TO DO:
// - Move the non-jump / non-treemove functions to Functions.cpp
// - General clean-up
//    - remove instances of tree.GetParameters() and declare a variable as such.
// - Update .h file


#include "TreeModifications.h"


//                           //
// Non-reversible tree moves // 
//                           //


bool Grow(Tree& tree, Cutpoints& cutpoints, Data& data, TreePrior& tree_prior, 
          const double& sigma, std::vector<size_t>& variable_inclusion_count, 
          std::vector<double>& variable_inclusion_prob, 
          ScaleMixture& scale_mixture, Random& random) {

  // Calculate the probability of growing a tree at a good bottom node
  std::vector<Tree*> leafs;  // Nodes that can potentially be split
  double PBx = GrowProbability(tree, cutpoints, tree_prior, leafs); // Probability of growing

  // Select a random bottom node from the list of good bottom nodes
  size_t g_node_index = floor(random.uniform() * leafs.size());
  Tree* g_node = leafs[g_node_index]; // The chosen node for the growth operation

  // Initialize variable for the cutpoint index
  size_t c = 0;

  // Retrieve the list of variables that the chosen node can split on
  std::vector<size_t> cut_variables;
  GetSplittableVariables(*g_node, cutpoints, cut_variables);

  // Select a variable to split on using inclusion probabilities
  random.SetInclusionWeights(variable_inclusion_prob); // Update weights based on probabilities
  size_t v = random.discrete(); // Draw the variable index to split on

  // Determine the range of cutpoints for the selected variable
  int lower_bound_index = 0;
  int upper_bound_index = cutpoints.values[v].size() - 1;

  // Check if the selected variable can actually be split
  if (!std::binary_search(cut_variables.begin(), cut_variables.end(), v)) {
    // If the selected variable cannot be split, use the cutpoint of the closest ancestor
    c = g_node->FindSameCut(v);
  } else {
    // If the variable is good for splitting, select a random cutpoint within the range
    c = lower_bound_index + 
        floor(random.uniform() * (upper_bound_index - lower_bound_index + 1));
  }

  // Calculate necessary probabilities for the Metropolis-Hastings acceptance ratio
  double Pbotx = 1.0 / leafs.size(); // Probability of choosing the node for splitting
  size_t dnx = g_node->NodeDepth();
  double PGnx = tree_prior.base / pow(1.0 + dnx, tree_prior.power); // Prior probability for growth

  double PGly, PGry; // Prior probabilities for the left and right child nodes
  if (cut_variables.size() > 1) {
    // If multiple variables can be split, assign prior probabilities
    PGly = tree_prior.base / pow(1.0 + dnx + 1.0, tree_prior.power); // Left child
    PGry = PGly; // Same for right child
  } else {
    // If only one variable is good for splitting, adjust probabilities based on cutpoint range
    PGly = ((int)(c - 1) < lower_bound_index) ? 
           0.0 : tree_prior.base / pow(1.0 + dnx + 1.0, tree_prior.power);
    PGry = (upper_bound_index < (int)(c + 1)) ? 
           0.0 : tree_prior.base / pow(1.0 + dnx + 1.0, tree_prior.power);
  }

  // Calculate the overall probability of the proposed operation
  double PDy = (leafs.size() > 1 || (PGry != 0 || PGly != 0)) ? 
               1.0 - tree_prior.p_GROW : 1.0;
  double Pnogy = (g_node->GetParent() == nullptr) ? 1.0 : 
                 (g_node->GetParent()->IsNog() ? 1.0 / tree.NumberOfNogs() : 
                  1.0 / (tree.NumberOfNogs() + 1.0));

  double pr = (PGnx * (1.0 - PGly) * (1.0 - PGry) * PDy * Pnogy) / 
              ((1.0 - PGnx) * Pbotx * PBx);

  // Compute sufficient statistics for the proposed child nodes
  size_t left_count, right_count; // Counts of observations in proposed nodes
  double left_residual, right_residual; // Sums of residuals in proposed nodes
  SufficientStatistics(tree, g_node, v, c, cutpoints, data, left_count, 
                       left_residual, right_count, right_residual);

  // Calculate the acceptance ratio (alpha) for the Metropolis-Hastings step
  double alpha = 0.0, log_alpha = 0.0;
  if ((left_count >= 5) && (right_count >= 5)) {  // Minimum node size requirement
    double lhl = LogPostLikelihood(left_count, left_residual, sigma, tree_prior.eta);
    double lhr = LogPostLikelihood(right_count, right_residual, sigma, tree_prior.eta);
    double lht = LogPostLikelihood(left_count + right_count, 
                                   left_residual + right_residual, 
                                   sigma, tree_prior.eta);
    log_alpha = log(pr) + (lhl + lhr - lht) + log(sigma);
    log_alpha = std::min(0.0, log_alpha); // Ensure log_alpha does not exceed 0
    alpha = exp(log_alpha);
  }

  // Metropolis-Hastings acceptance step
  Parameters par_left(g_node->GetParameters()), par_right(g_node->GetParameters());
  if (alpha > 0 && log(random.uniform()) < log_alpha) {

    // If accepted, draw new mu values for the left and right nodes
    par_left.SetParameters(0, 
        DrawMuOneLeave(left_count, left_residual, tree_prior.eta, sigma, random));
    par_right.SetParameters(0, 
        DrawMuOneLeave(right_count, right_residual, tree_prior.eta, sigma, random));
    
    // Grow the tree with the new child nodes
    tree.GrowChildren(g_node, v, c, par_left, par_right);
    variable_inclusion_count[v]++;
    return true;  // Grow step accepted
  } else {
    return false; // Grow step rejected
  }
}


bool Prune(Tree& tree, Cutpoints& cutpoints, Data& data, TreePrior& tree_prior, 
           const double& sigma, std::vector<size_t>& variable_inclusion_count, 
           std::vector<double>& variable_inclusion_prob, 
           ScaleMixture& scale_mixture, Random& random) {

  // Calculate the probability of growing a tree at a good bottom node
  std::vector<Tree*> leafs;  // Nodes available for splitting
  double PBx = GrowProbability(tree, cutpoints, tree_prior, leafs);

  // Collect all NOG nodes (non-terminal nodes with both children being terminal)
  std::vector<Tree*> nog_nodes;
  tree.CollectNogs(nog_nodes);

  // Randomly select one NOG node for pruning
  size_t p_node_id = floor(random.uniform() * nog_nodes.size());
  Tree* p_node = nog_nodes[p_node_id];

  // Calculate probabilities needed for the Metropolis-Hastings acceptance ratio

  // Probability that the selected NOG node would grow
  double PGny = tree_prior.base / pow(1.0 + p_node->NodeDepth(), 
                                       tree_prior.power);

  // Probabilities of growth for the left and right child nodes
  double PGlx = ProbNodeGrows(*(p_node->GetLeft()), cutpoints, tree_prior);
  double PGrx = ProbNodeGrows(*(p_node->GetRight()), cutpoints, tree_prior);

  // Probability of proposing a birth move at the selected node
  double PBy = (p_node->GetParent() == nullptr) ? 1.0 : tree_prior.p_GROW;

  // Probability of choosing the selected NOG node as a bottom node to split on
  int ngood = nog_nodes.size();
  if (Splittable(*(p_node->GetLeft()), cutpoints)) --ngood;
  if (Splittable(*(p_node->GetRight()), cutpoints)) --ngood;
  ++ngood;  // Include the selected node itself as splittable
  double Pboty = 1.0 / ngood;

  // Probability of a death move and selecting the chosen NOG node
  double PDx = 1.0 - PBx; // Probability of death move
  double Pnogx = 1.0 / nog_nodes.size(); // Probability of choosing NOG node

  // Calculate part of the Metropolis ratio from proposal and prior
  double pr = ((1.0 - PGny) * PBy * Pboty) / 
              (PGny * (1.0 - PGlx) * (1.0 - PGrx) * PDx * Pnogx);

  // Compute sufficient statistics for the left and right child nodes
  size_t count_left, count_right;
  double residual_left, residual_right;
  SufficientStatistics(tree, p_node->GetLeft(), p_node->GetRight(), 
                       cutpoints, data, count_left, residual_left, 
                       count_right, residual_right);

  // Calculate the log-likelihood ratio
  double lhl = LogPostLikelihood(count_left, residual_left, sigma, 
                                  tree_prior.eta);
  double lhr = LogPostLikelihood(count_right, residual_right, sigma, 
                                  tree_prior.eta);
  double lht = LogPostLikelihood(count_left + count_right, 
                                  residual_left + residual_right, 
                                  sigma, tree_prior.eta);
  double log_likelihood_ratio = lht - lhl - lhr - log(sigma);

  // Calculate the acceptance ratio (alpha)
  double log_alpha = log(pr) + log_likelihood_ratio;
  log_alpha = std::min(0.0, log_alpha); // Ensure log_alpha does not exceed 0

  // Perform the Metropolis-Hastings acceptance step
  if (log(random.uniform()) < log_alpha) {
    Parameters par(tree.GetParameters());
    par.SetParameters(0, DrawMuOneLeave(count_left + count_right, 
                                         residual_left + residual_right, 
                                         tree_prior.eta, sigma, random));
    variable_inclusion_count[p_node->GetSplitVar()]--;
    tree.KillChildren(p_node, par);
    return true;  // Pruning accepted
  } else {
    return false; // Pruning rejected
  }
}


bool Change(Tree& tree, Cutpoints& cutpoints, Data& data, TreePrior& tree_prior, 
            const double& sigma, std::vector<size_t>& variable_inclusion_count, 
            std::vector<double>& variable_inclusion_prob, 
            ScaleMixture& scale_mixture, Random& random) {

  // Collect all NOG nodes (non-terminal nodes) and randomly choose one to change
  std::vector<Tree*> nog_nodes;
  tree.CollectNogs(nog_nodes);
  size_t c_node_id = floor(random.uniform() * nog_nodes.size());
  Tree* c_node = nog_nodes[c_node_id]; // Selected NOG node for modification

  // Compute sufficient statistics for the current split
  size_t count_left_current, count_right_current; // Counts in the child nodes
  double residual_left_current, residual_right_current; // Sums at child nodes
  SufficientStatistics(tree, c_node->GetLeft(), c_node->GetRight(), 
                       cutpoints, data, count_left_current, 
                       residual_left_current, count_right_current, 
                       residual_right_current);

  // Compute log-likelihood ratios for the current state
  double log_lik_left_current = LogPostLikelihood(count_left_current, 
                                                   residual_left_current, 
                                                   sigma, tree_prior.eta);
  double log_lik_right_current = LogPostLikelihood(count_right_current, 
                                                    residual_right_current, 
                                                    sigma, tree_prior.eta);

  // Get all variables that can be split and draw a new variable for the split
  std::vector<size_t> variables; 
  GetSplittableVariables(*c_node, cutpoints, variables);
  random.SetInclusionWeights(variable_inclusion_prob);    
  size_t proposed_split_var = random.discrete(); // Draw new split variable

  // Draw a cut value for the proposed splitting variable
  size_t proposed_cut_val;
  int lower_bound_index = 0;
  int upper_bound_index = cutpoints.values[proposed_split_var].size() - 1;

  // Check if the selected variable can be split
  if (!std::binary_search(variables.begin(), variables.end(), proposed_split_var)) {
    // If the variable cannot be split, use the cutpoint of the closest ancestor
    proposed_cut_val = c_node->FindSameCut(proposed_split_var);
  } else {
    // If the variable is valid for splitting, select a random cutpoint
    c_node->PossibleCuts(proposed_split_var, &lower_bound_index, 
                         &upper_bound_index);
    proposed_cut_val = lower_bound_index + 
                       floor(random.uniform() * (upper_bound_index - lower_bound_index + 1));
  }

  // Compute sufficient statistics for the proposed split
  size_t count_left_proposed, count_right_proposed; // Counts for proposed nodes
  double residual_left_proposed, residual_right_proposed; // Sums for proposed nodes
  SufficientStatistics(tree, c_node, proposed_split_var, proposed_cut_val, 
                       cutpoints, data, count_left_proposed, 
                       residual_left_proposed, count_right_proposed, 
                       residual_right_proposed);

  // Compute log-likelihood ratios for the proposed state
  double log_lik_left_proposed = LogPostLikelihood(count_left_proposed, 
                                                    residual_left_proposed, 
                                                    sigma, tree_prior.eta);
  double log_lik_right_proposed = LogPostLikelihood(count_right_proposed, 
                                                     residual_right_proposed, 
                                                     sigma, tree_prior.eta);

  // Calculate the acceptance ratio (alpha) for the Metropolis-Hastings step
  double log_alpha = log_lik_left_proposed + log_lik_right_proposed - 
                     (log_lik_left_current + log_lik_right_current);
  log_alpha = std::min(0.0, log_alpha); // Ensure log_alpha does not exceed 0

  // Perform the Metropolis-Hastings acceptance step
  if (log(random.uniform()) < log_alpha) {
    variable_inclusion_count[c_node->GetSplitVar()]--;
    variable_inclusion_count[proposed_split_var]++;
    c_node->SetSplitVar(proposed_split_var);
    c_node->SetCutVal(proposed_cut_val);
    return true; // Change accepted
  } else {
    return false; // Change rejected
  }
}


//                            //
// Reversible Jump tree moves // 
//                            //


bool RJ_Grow(Tree& tree, Cutpoints& cutpoints, Data& data, TreePrior& tree_prior, 
             const double& sigma, const double& omega, std::vector<size_t>& variable_inclusion_count, 
             std::vector<double>& variable_inclusion_prob, const size_t& delayed_proposal,
             ScaleMixture& scale_mixture, Random& random) {

  // Collect all leaf nodes to potentially grow from
  std::vector<Tree*> leafs;  
  tree.CollectLeafs(leafs);
  size_t number_of_leafs = leafs.size();
  
  // Select a random bottom node from the list of leaf nodes
  size_t g_node_index = floor(random.uniform() * number_of_leafs);
  Tree* g_node = leafs[g_node_index]; // The selected leaf node for the growth

  // Initialize the proposed cutpoint variable
  size_t c = 0;

  // Retrieve the variables that the selected node can split on
  std::vector<size_t> cut_variables;
  GetSplittableVariables(*g_node, cutpoints, cut_variables);

  // Use inclusion probabilities to select a variable for splitting
  random.SetInclusionWeights(variable_inclusion_prob);
  size_t variable_index = floor(random.uniform() * cut_variables.size());
  size_t v = cut_variables[variable_index]; // The proposed variable to split on

  // Determine the range of cutpoints for the selected variable
  int lower_bound_index = 0;
  int upper_bound_index = cutpoints.values[v].size() - 1;

  // Draw a cutpoint for the proposed variable
  c = lower_bound_index + 
      floor(random.uniform() * (upper_bound_index - lower_bound_index + 1));

  // Compute sufficient statistics for the proposed split
  size_t left_count, right_count;  // Counts for the left and right child nodes
  double left_residual, right_residual; // Sums of outcomes for proposed nodes
  SufficientStatistics(tree, g_node, v, c, cutpoints, data, left_count, 
                       left_residual, right_count, right_residual);

  // Initialize parameters for the proposed child nodes
  Parameters par_left(g_node->GetParameters()), 
             par_right(g_node->GetParameters());

  // Calculate the acceptance ratio (alpha) based on sufficient statistics
  double alpha = 0.0, log_alpha = 0.0;
  if ((left_count >= 5) && (right_count >= 5)) {  // Ensure minimum node size
    size_t depth = g_node->NodeDepth(); // Depth of the selected node
    double grow_prob = (g_node->GetParent() == nullptr) ? 
                       1.0 : tree_prior.p_GROW;

    size_t number_of_nogs = (g_node->GetParent() == nullptr) ? 
                             1.0 : (g_node->GetParent()->IsNog() ? 
                             tree.NumberOfNogs() : (tree.NumberOfNogs() + 1.0));

    // Propose new mu values for the left and right nodes
    for (size_t p = 0; p < delayed_proposal; p++) {
      scale_mixture.Propose(par_left, tree.GetParameters(), left_residual, 
                            left_count, sigma, omega, random);
      scale_mixture.Propose(par_right, tree.GetParameters(), right_residual, 
                            right_count, sigma, omega, random);
    }

    // Compute the log-likelihood ratio for the proposed split
    double log_likelihood_ratio = scale_mixture.LogLikelihood(par_left, 
        tree.GetParameters(), left_residual, left_count, sigma) + 
        scale_mixture.LogLikelihood(par_right, tree.GetParameters(), 
        right_residual, right_count, sigma) - 
        scale_mixture.LogLikelihood(g_node->GetParameters(), 
        tree.GetParameters(), left_residual + right_residual, 
        left_count + right_count, sigma); 

    // Calculate the log move ratio and tree ratio for GROW move
    double log_move_ratio = LogMoveRatio(number_of_nogs, number_of_leafs, 
                                          grow_prob, tree_prior.p_PRUNE);
    double log_tree_ratio = LogTreeRatio_GROW(depth, tree_prior);
    
    // Calculate the prior ratio based on the new parameters
    double log_prior_ratio = scale_mixture.LogPrior(par_left, 
        tree.GetParameters(), sigma, omega) + scale_mixture.LogPrior(par_right, 
        tree.GetParameters(), sigma, omega) - 
        scale_mixture.LogPrior(g_node->GetParameters(), 
        tree.GetParameters(), sigma, omega);

    // Compute the log of the proposal density ratio
    double log_proposition_ratio = scale_mixture.LogProposeDensity(
        g_node->GetParameters(), tree.GetParameters(), left_residual + 
        right_residual, left_count + right_count, sigma, omega) - 
        scale_mixture.LogProposeDensity(par_left, tree.GetParameters(), 
        left_residual, left_count, sigma, omega) - 
        scale_mixture.LogProposeDensity(par_right, tree.GetParameters(), 
        right_residual, right_count, sigma, omega);

    // Combine all components to compute log_alpha
    log_alpha = log_likelihood_ratio + log_move_ratio + log_tree_ratio + 
                log_proposition_ratio + log_prior_ratio;
    log_alpha = std::min(0.0, log_alpha); // Ensure log_alpha does not exceed 0
    alpha = exp(log_alpha); // Calculate acceptance probability
  }

  // Perform the Metropolis-Hastings acceptance step
  if (alpha > 0 && log(random.uniform()) < log_alpha) {
    // If accepted, update the tree with new parameters
    tree.GrowChildren(g_node, v, c, par_left, par_right);
    variable_inclusion_count[v]++;
    return true; // Grow step accepted
  } else {
    return false; // Grow step rejected
  }
}


bool RJ_Prune(Tree& tree, Cutpoints& cutpoints, Data& data, TreePrior& tree_prior, 
              const double& sigma, const double& omega, std::vector<size_t>& variable_inclusion_count, 
              std::vector<double>& variable_inclusion_prob, const size_t& delayed_proposal,
              ScaleMixture& scale_mixture, Random& random) {

  // Collect all prunable nodes, which are the NOG nodes
  std::vector<Tree*> nog_nodes;
  tree.CollectNogs(nog_nodes);

  // Randomly select one NOG node for pruning
  size_t p_node_id = floor(random.uniform() * nog_nodes.size());
  Tree* p_node = nog_nodes[p_node_id];

  // Compute sufficient statistics for the current child nodes of the selected node
  size_t left_count, right_count; // Counts for left and right child nodes
  double left_residual, right_residual; // Sums of outcomes for child nodes
  SufficientStatistics(tree, p_node->GetLeft(), p_node->GetRight(), cutpoints, 
                       data, left_count, left_residual, right_count, 
                       right_residual);

  // Retrieve global parameters from the tree
  Parameters global_parameters = tree.GetParameters();

  // Compute log-likelihood for the current split
  double log_lik_current = scale_mixture.LogLikelihood(p_node->GetLeft()->GetParameters(), 
        global_parameters, left_residual, left_count, sigma) + 
        scale_mixture.LogLikelihood(p_node->GetRight()->GetParameters(), 
        global_parameters, right_residual, right_count, sigma);

  // Compute log prior density at current parameter values
  double prior_current = scale_mixture.LogPrior(p_node->GetLeft()->GetParameters(), 
        global_parameters, sigma, omega) + scale_mixture.LogPrior(p_node->GetRight()->GetParameters(), 
        global_parameters, sigma, omega);

  // Compute log density at current parameter values
  double proposed_density_current = scale_mixture.LogProposeDensity(p_node->GetLeft()->GetParameters(), 
        global_parameters, left_residual, left_count, sigma, omega) + 
        scale_mixture.LogProposeDensity(p_node->GetRight()->GetParameters(), 
        global_parameters, right_residual, right_count, sigma, omega);

  // Propose a new value for the local parameters of the prunable node

  for (size_t p = 0; p < delayed_proposal; p++) {
    scale_mixture.Propose(p_node->GetParameters(), global_parameters, 
                          left_residual + right_residual, 
                          left_count + right_count, sigma, omega, random);
  }

  // Compute log-likelihood for the proposed split
  double log_lik_new = scale_mixture.LogLikelihood(p_node->GetParameters(), 
        global_parameters, left_residual + right_residual, 
        left_count + right_count, sigma);

  // Compute log prior density at proposed parameter values
  double prior_new = scale_mixture.LogPrior(p_node->GetParameters(), global_parameters, sigma, omega);

  // Compute proposed density at current parameter values
  double proposed_density_new = scale_mixture.LogProposeDensity(p_node->GetParameters(), 
        global_parameters, left_residual + right_residual, 
        left_count + right_count, sigma, omega);

  // Compute the probability that the selected node will grow
  double grow_prob = (p_node == nullptr) ? 1.0 : tree_prior.p_GROW;

  // Compute the move ratio and tree ratio for the pruning operation
  size_t number_of_nogs = tree.NumberOfNogs();
  size_t number_of_leafs = tree.NumberOfLeafs();

  double log_move_ratio = -LogMoveRatio(number_of_nogs, number_of_leafs, 
                                         grow_prob, tree_prior.p_PRUNE);
  double log_tree_ratio = LogTreeRatio_PRUNE(p_node, tree_prior, cutpoints);

  // Calculate the acceptance ratio (alpha)
  double log_alpha = log_move_ratio + log_tree_ratio + log_lik_new + 
                     prior_new - proposed_density_new - 
                     (log_lik_current + prior_current - proposed_density_current);
  log_alpha = std::min(0.0, log_alpha); // Ensure log_alpha does not exceed 0

  // Perform the Metropolis-Hastings acceptance step
  if (log(random.uniform()) < log_alpha) {
    variable_inclusion_count[p_node->GetSplitVar()]--;
    tree.KillChildren(p_node); // Remove the children of the pruned node
    return true;  // Pruning accepted
  } else {
    return false; // Pruning rejected
  }
}


bool RJ_Change(Tree& tree, Cutpoints& cutpoints, Data& data, TreePrior& tree_prior, 
               const double& sigma, const double& omega, std::vector<size_t>& variable_inclusion_count, 
               std::vector<double>& variable_inclusion_prob, const size_t& delayed_proposal,
               ScaleMixture& scale_mixture, Random& random) {

  // Collect all NOG nodes (non-terminal nodes)
  std::vector<Tree*> nog_nodes;
  tree.CollectNogs(nog_nodes);

  // Randomly select one NOG node to change its split variable and cut value
  size_t c_node_id = floor(random.uniform() * nog_nodes.size());
  Tree* c_node = nog_nodes[c_node_id];

  // Compute sufficient statistics for the current child nodes of the selected node
  size_t count_left_current, count_right_current; // Counts in left and right child nodes
  double residual_left_current, residual_right_current; // Sums at child nodes
  SufficientStatistics(tree, c_node->GetLeft(), c_node->GetRight(), 
                       cutpoints, data, count_left_current, 
                       residual_left_current, count_right_current, 
                       residual_right_current);

  // Retrieve global parameters from the tree
  Parameters global_parameters = tree.GetParameters();

  // Compute log-likelihood for the current split
  double log_lik_current = scale_mixture.LogLikelihood(c_node->GetLeft()->GetParameters(), 
        global_parameters, residual_left_current, count_left_current, sigma) + 
        scale_mixture.LogLikelihood(c_node->GetRight()->GetParameters(), 
        global_parameters, residual_right_current, count_right_current, sigma);

  // Compute log prior density at current parameter values
  double prior_current = scale_mixture.LogPrior(c_node->GetLeft()->GetParameters(), 
        global_parameters, sigma, omega) + scale_mixture.LogPrior(c_node->GetRight()->GetParameters(), 
        global_parameters, sigma, omega);

  // Compute log density at current parameter values
  double proposed_density_current = scale_mixture.LogProposeDensity(c_node->GetLeft()->GetParameters(), 
        global_parameters, residual_left_current, count_left_current, sigma, omega) + 
        scale_mixture.LogProposeDensity(c_node->GetRight()->GetParameters(), 
        global_parameters, residual_right_current, count_right_current, sigma, omega);

  // Select all variables that can be split and draw a new variable to split on
  std::vector<size_t> variables; 
  GetSplittableVariables(*c_node, cutpoints, variables);
  random.SetInclusionWeights(variable_inclusion_prob); // Update weights based on inclusion probabilities
  size_t proposed_split_var = random.discrete(); // Draw new split variable

  // Draw a cut value for the proposed splitting variable
  size_t proposed_cut_val;
  int lower_bound_index = 0;
  int upper_bound_index = cutpoints.values[proposed_split_var].size() - 1;

  // Check if the selected variable can be split
  if (!std::binary_search(variables.begin(), variables.end(), proposed_split_var)) {
    // If the variable cannot be split, use the cutpoint of the closest ancestor
    proposed_cut_val = c_node->FindSameCut(proposed_split_var);
  } else {
    // If the variable is valid for splitting, select a random cutpoint
    c_node->PossibleCuts(proposed_split_var, &lower_bound_index, &upper_bound_index);
    proposed_cut_val = lower_bound_index + 
                       floor(random.uniform() * (upper_bound_index - lower_bound_index + 1));
  }

  // Compute sufficient statistics for the proposed split
  size_t count_left_proposed, count_right_proposed; // Counts for proposed child nodes
  double residual_left_proposed, residual_right_proposed; // Sums for proposed child nodes
  SufficientStatistics(tree, c_node, proposed_split_var, proposed_cut_val, 
                       cutpoints, data, count_left_proposed, 
                       residual_left_proposed, count_right_proposed, 
                       residual_right_proposed);

  // Initialize parameters for the proposed child nodes
  Parameters par_left(global_parameters), par_right(global_parameters);
  for (size_t p = 0; p < delayed_proposal; p++) {
    scale_mixture.Propose(par_left, global_parameters, residual_left_proposed, 
                          count_left_proposed, sigma, omega, random);
    scale_mixture.Propose(par_right, global_parameters, residual_right_proposed, 
                          count_right_proposed, sigma, omega, random);
  }

  // Compute log-likelihood for the proposed split
  double log_lik_new = scale_mixture.LogLikelihood(par_left, global_parameters, 
        residual_left_proposed, count_left_proposed, sigma) + 
        scale_mixture.LogLikelihood(par_right, global_parameters, 
        residual_right_proposed, count_right_proposed, sigma);

  // Compute log prior density at proposed parameter values
  double prior_new = scale_mixture.LogPrior(par_left, global_parameters, sigma, omega) + 
                     scale_mixture.LogPrior(par_right, global_parameters, sigma, omega);

  // Compute proposed density at current parameter values
  double proposed_density_new = scale_mixture.LogProposeDensity(par_left, 
        global_parameters, residual_left_proposed, count_left_proposed, sigma, omega) + 
        scale_mixture.LogProposeDensity(par_right, global_parameters, 
        residual_right_proposed, count_right_proposed, sigma, omega);

  // Calculate the acceptance ratio (alpha) for the Metropolis-Hastings step
  double log_alpha = log_lik_new + prior_new - proposed_density_new - 
                     (log_lik_current + prior_current - proposed_density_current);
  log_alpha = std::min(0.0, log_alpha); // Ensure log_alpha does not exceed 0

  // Perform the Metropolis-Hastings acceptance step
  if (log(random.uniform()) < log_alpha) {
    // If accepted, update the tree with new parameters
    variable_inclusion_count[c_node->GetSplitVar()]--;
    variable_inclusion_count[proposed_split_var]++;
    c_node->SetSplitVar(proposed_split_var);
    c_node->SetCutVal(proposed_cut_val);
    c_node->GetLeft()->SetParameters(par_left); // Update left child parameters
    c_node->GetRight()->SetParameters(par_right); // Update right child parameters
    return true; // Change accepted
  } else {
    return false; // Change rejected
  }
}


//                                            //
// Wrapper of tree moves and reversible jumps //
//                                            //


bool ReversibleJump(Tree& tree, Cutpoints& cutpoints, Data& data, 
                    TreePrior& tree_prior, const double& sigma, const double& omega, 
                    std::vector<size_t>& nodes, 
                    std::vector<double>& variable_inclusion_prob, 
                    const size_t& delayed_proposal,
                    bool reversible_jump, ScaleMixture& scale_mixture, 
                    Random& random) {

  // Set the probabilities for each action: grow, prune, change
  double prob_GROW = tree_prior.p_GROW;
  double prob_PRUNE = tree_prior.p_PRUNE;
  double U;

  // Indicator for whether to use Reversible Jump MCMC or standard BART
  bool RJMCMC = reversible_jump;

  // If using Reversible Jump MCMC
  if (RJMCMC) {
    // For a stump (a tree with only one node), we can only grow
    if (tree.TreeSize() == 1) {
      return RJ_Grow(tree, cutpoints, data, tree_prior, sigma, omega, nodes, 
                     variable_inclusion_prob, delayed_proposal, scale_mixture, random);
    } else {
      // Randomly select one of the actions based on the specified probabilities
      U = random.uniform();

      // Attempt to grow the tree
      if (U < prob_GROW) {
        return RJ_Grow(tree, cutpoints, data, tree_prior, sigma, omega, nodes, 
                       variable_inclusion_prob, delayed_proposal, scale_mixture, random);
      } 
      // Attempt to prune the tree
      else if (U < prob_GROW + prob_PRUNE) {
        return RJ_Prune(tree, cutpoints, data, tree_prior, sigma, omega, nodes, 
                        variable_inclusion_prob, delayed_proposal, scale_mixture, random);
      } 
      // Attempt to change the tree
      else {
        return RJ_Change(tree, cutpoints, data, tree_prior, sigma, omega, nodes, 
                         variable_inclusion_prob, delayed_proposal, scale_mixture, random);
      }
    }
  } 
  // If not using Reversible Jump MCMC
  else {
    // For a stump, we can only grow
    if (tree.TreeSize() == 1) {
      return Grow(tree, cutpoints, data, tree_prior, sigma, nodes, 
                  variable_inclusion_prob, scale_mixture, random);
    } else {
      // Randomly select one of the actions based on the specified probabilities
      U = random.uniform();

      // Attempt to grow the tree
      if (U < prob_GROW) {
        return Grow(tree, cutpoints, data, tree_prior, sigma, nodes, 
                    variable_inclusion_prob, scale_mixture, random);
      } 
      // Attempt to prune the tree
      else if (U < prob_GROW + prob_PRUNE) {
        return Prune(tree, cutpoints, data, tree_prior, sigma, nodes, 
                     variable_inclusion_prob, scale_mixture, random);
      } 
      // Attempt to change the tree
      else {
        return Change(tree, cutpoints, data, tree_prior, sigma, nodes, 
                      variable_inclusion_prob, scale_mixture, random);
      }
    }
  }
}
