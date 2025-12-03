#ifndef GUARD_StanTreeFunctions_h
#define GUARD_StanTreeFunctions_h

#include "StanTree.h"


//fit StanTree at matrix of x, matrix is stacked columns x[i,j] is *(x+p*i+j)
void fit(StanTree& t, xinfo& xi, size_t p, size_t n, double *x,  double* fv);

//does a (bottom) node have variables you can split on?
bool cansplit(StanTree::StanTree_p n, xinfo& xi);

//find variables n can split on, put their indices in goodvars
void getgoodvars(StanTree::StanTree_p n, xinfo& xi,  std::vector<size_t>& goodvars);

#endif
