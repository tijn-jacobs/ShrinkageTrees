#ifndef GUARD_StanForestFunctions_h
#define GUARD_StanForestFunctions_h

#include "StanTree.h"
#include "StanTreeFunctions.h"
#include "Info.h"
#include <algorithm>


//make xinfo = cutpoints
void makexinfo(size_t p, size_t n, double *x, xinfo& xi, int* nc);

//compute prob of a birth, goodbots will contain all the good bottom nodes
double getpb(StanTree& t, xinfo& xi, pinfo& pi, StanTree::npv& goodbots);

//compute n and \sum y_i for left and right give bot and v,c
void getsuff(StanTree& x, StanTree::StanTree_p nx, size_t v, size_t c, xinfo& xi, dinfo& di, size_t& nl, double& syl, size_t& nr, double& syr);

//lh, replacement for lil that only depends on sum y.
double lh(size_t n, double sy, double sigma, double eta);

//get prob a node grows, 0 if no good vars, else a/(1+d)^b
double pgrow(StanTree::StanTree_p n, xinfo& xi, pinfo& pi);

//compute n and \sum y_i for left and right bots
void getsuff(StanTree& x, StanTree::StanTree_p l, StanTree::StanTree_p r, xinfo& xi, dinfo& di, size_t& nl, double& syl, size_t& nr, double& syr);

//get sufficients stats for all bottom nodes, this way just loop through all the data once.
void allsuff(StanTree& x, xinfo& xi, dinfo& di, StanTree::npv& leaves, std::vector<size_t>& nv, std::vector<double>& syv);

// draw all the bottom node mu's
void drmu(StanTree& t, xinfo& xi, dinfo& di, pinfo& pi, double sigma, Random& random);

//birth proposal
void bprop(StanTree& x, xinfo& xi, pinfo& pi, StanTree::npv& goodbots, double& PBx, StanTree::StanTree_p& nx, size_t& v, size_t& c, double& pr, std::vector<size_t>& nv, std::vector<double>& pv, bool aug, Random& random);

// death proposal
void dprop(StanTree& x, xinfo& xi, pinfo& pi, StanTree::npv& goodbots, double& PBx, StanTree::StanTree_p& nx, double& pr, Random& random);

//draw one mu from post 
double drawnodemu(size_t n, double sy, double eta, double sigma, Random& random);

//draw variable splitting probabilities from Dirichlet (Linero, 2018)
void draw_s(std::vector<size_t>& nv, std::vector<double>& lpv, double& theta, Random& random);

//draw Dirichlet sparsity parameter from posterior using grid
void draw_theta0(bool const_theta, double& theta, std::vector<double>& lpv,
		 double a, double b, double rho, Random& random);
#endif
