#ifndef GUARD_Dirichlet_h
#define GUARD_Dirichlet_h

#include "Forest.h"

#include <vector>
using std::vector; // Double?!?!?!

double LogSumExp(const vector<double>& values);

vector<double> LogDirichlet(const vector<double>& dirichlet_shape,
                            Random& random);

void DrawSplitProbs(vector<size_t>& variable_inclusion_count,
                    vector<double>& log_split_probs,
                    double& alpha_dirichlet,
                    Random& random);

void DrawDirichletAlpha(bool const_alpha,
                        double& alpha_dirichlet,
                        const vector<double>& log_split_probs,
                        vector<double>& variable_inclusion_prob,
                        double a_dirichlet,
                        double b_dirichlet,
                        double rho_dirichlet,
                        Random& random);

#endif // GUARD_Dirichlet_h
