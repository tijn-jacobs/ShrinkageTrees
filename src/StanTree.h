#ifndef GUARD_StanTree_h
#define GUARD_StanTree_h

#include "Prerequisites.h"


//xinfo xi, then xi[v][c] is the c^{th} cutpoint for variable v.
//left if x[v] < xi[v][c]
typedef std::vector<double> vec_d; //double vector
typedef std::vector<vec_d> xinfo; //vector of vectors, will be split rules


//info contained in a node, used by input operator
struct node_info {
   std::size_t id; // node id
   std::size_t split_var;  // variable
   std::size_t cut_val;  // cut point
   double step_height;   // step_height
};


class StanTree {
public:
   //friends--------------------
   friend std::istream& operator>>(std::istream&, StanTree&);
   //typedefs--------------------
   typedef StanTree* StanTree_p;
   typedef const StanTree* StanTree_cp;
   typedef std::vector<StanTree_p> npv; 
   typedef std::vector<StanTree_cp> cnpv;
   //contructors,destructors--------------------
   StanTree(): step_height(0.0),split_var(0),cut_val(0),p(0),l(0),r(0) {}
   StanTree(const StanTree& n): step_height(0.0),split_var(0),cut_val(0),p(0),l(0),r(0) {cp(this,&n);}
   StanTree(double istep_height): step_height(istep_height),split_var(0),cut_val(0),p(0),l(0),r(0) {}
   void tonull(); //like a "clear", null StanTree has just one node
   ~StanTree() {tonull();}
   //operators----------
   StanTree& operator=(const StanTree&);
   //interface--------------------
   //set
   void settheta(double step_height) {this->step_height=step_height;}
   void setv(size_t split_var) {this->split_var = split_var;}
   void setc(size_t cut_val) {this->cut_val = cut_val;}
   //get
   double gettheta() const {return step_height;}
   size_t GetSplitVar() const {return split_var;}
   size_t GetCutVal() const {return cut_val;}
   StanTree_p getp() {return p;}  
   StanTree_p getl() {return l;}
   StanTree_p getr() {return r;}
   //StanTree functions--------------------
   StanTree_p getptr(size_t nid); //get node pointer from node id, 0 if not there
   void pr(bool pc=true); //to screen, pc is "print children"
   size_t StanTreesize(); //number of nodes in StanTree
   size_t nnogs();    //number of nog nodes (no grandchildren nodes)
   size_t nbots();    //number of bottom nodes
   bool birth(size_t nid, size_t split_var, size_t cut_val, double step_heightl, double step_heightr);
   bool death(size_t nid, double step_height);
   void birthp(StanTree_p np, size_t split_var, size_t cut_val, double step_heightl, double step_heightr);
   void deathp(StanTree_p nb, double step_height);
   void getbots(npv& bv);         //get bottom nodes
   void getnogs(npv& nv);         //get nog nodes (no granchildren)
   void getnodes(npv& v);         //get vector of all nodes
   void getnodes(cnpv& v) const;  //get vector of all nodes (const)
   StanTree_p bn(double *x,xinfo& xi); //find Bottom Node
   void rg(size_t split_var, int* L, int* U); //recursively find region [L,U] for var v
   //node functions--------------------
   size_t nid() const; //nid of a node
   size_t depth();  //depth of a node
   char ntype(); //node type t:top, b:bot, n:no grandchildren i:interior (t can be b)
   bool isnog();
   size_t getbadcut(size_t split_var);
#ifndef NoRcpp   
  Rcpp::List StanTree2list(xinfo& xi, double center=0., double scale=1.); // create an efficient list from a single StanTree
  //StanTree list2StanTree(Rcpp::List&, xinfo& xi); // create a StanTree from a list and an xinfo  
  Rcpp::IntegerVector StanTree2count(size_t nvar); // for one StanTree, count the number of branches for each variable
#endif
private:
   double step_height; //univariate double parameter
   //rule: left if x[cut_val] < xinfo[split_var][cut_val]
   size_t split_var;
   size_t cut_val;
   //StanTree structure
   StanTree_p p; //parent
   StanTree_p l; //left child
   StanTree_p r; //right child
   //utiity functions
   void cp(StanTree_p n,  StanTree_cp o); //copy StanTree
};
std::istream& operator>>(std::istream&, StanTree&);
std::ostream& operator<<(std::ostream&, const StanTree&);

#endif
