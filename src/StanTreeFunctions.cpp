#include "StanTreeFunctions.h"


//fit StanTree at matrix of x, matrix is stacked columns x[i,j] is *(x+p*i+j)
void fit(StanTree& t, xinfo& xi, size_t p, size_t n, double *x,  double* fv)
{
   StanTree::StanTree_p bn;
   for(size_t i=0;i<n;i++) {
      bn = t.bn(x+i*p,xi);
      fv[i] = bn->gettheta();
   }
}

//does this bottom node n have any variables it can split on.
bool cansplit(StanTree::StanTree_p n, xinfo& xi)
{
   int L,U;
   bool v_found = false; //have you found a variable you can split on
   size_t v=0;
   while(!v_found && (v < xi.size())) { //invar: splitvar not found, vars left
      L=0; U = xi[v].size()-1;
      n->rg(v,&L,&U);
      if(U>=L) v_found=true;
      v++;
   }
   return v_found;
}


//find variables n can split on, put their indices in goodvars
void getgoodvars(StanTree::StanTree_p n, xinfo& xi,  std::vector<size_t>& goodvars)
{
   goodvars.clear();
   int L,U;
   for(size_t v=0;v!=xi.size();v++) {//try each variable
      L=0; U = xi[v].size()-1;
      n->rg(v,&L,&U);
      if(U>=L) goodvars.push_back(v);
   }
}
