#ifndef GUARD_Dirichlet_h
#define GUARD_Dirichlet_h

#include "Forest.h"

#include <vector>
using std::vector; // Double?!?!?!

double LogSumExp(const vector<double>& values);

vector<double> LogDirichlet(const vector<double>& dirichlet_shape,
                            Random& random);

void DrawSplitProbs(vector<size_t>& variable_inclusion_count,
                    vector<double>& variable_inclusion_prob,
                    double& alpha_dirichlet,
                    Random& random);

                    /*
void DrawDirichletAlpha(bool const_alpha,
                        double& alpha_dirichlet,
                        vector<double>& variable_inclusion_prob,
                        double a_dirichlet,
                        double b_dirichlet,
                        double rho_dirichlet,
                        Random& random);
                        */

void DrawDirichletAlpha(bool const_theta, double& theta, std::vector<double>& lpv,
		 double a, double b, double rho, Random& random);

double log_sum_exp(std::vector<double>& v);
#endif // GUARD_Dirichlet_h
