
#include "StanTree.h"


// node id
size_t StanTree::nid() const 
{
   if(!p) return 1; //if you don't have a parent, you are the top
   if(this==p->l) return 2*(p->nid()); //if you are a left child
   else return 2*(p->nid())+1; //else you are a right child
}

StanTree::StanTree_p StanTree::getptr(size_t nid)
{
   if(this->nid() == nid) return this; //found it
   if(l==0) return 0; //no children, did not find it
   StanTree_p lp = l->getptr(nid);
   if(lp) return lp; //found on left
   StanTree_p rp = r->getptr(nid);
   if(rp) return rp; //found on right
   return 0; //never found it
}

//add children to  bot node nid
bool StanTree::birth(size_t nid, size_t split_var, size_t cut_val, double step_heightl, double step_heightr)
{
   StanTree_p np = getptr(nid);
   if(np==0) {
      cout << "error in birth: bottom node not found\n";
      return false; //did not find note with that nid
   }
   if(np->l!=0) {
      cout << "error in birth: found node has children\n";
      return false; //node is not a bottom node
   }

   //add children to bottom node np
   StanTree_p l = new StanTree;
   l->step_height=step_heightl;
   StanTree_p r = new StanTree;
   r->step_height=step_heightr;
   np->l=l;
   np->r=r;
   np->split_var = split_var; np->cut_val=cut_val;
   l->p = np;
   r->p = np;

   return true;
}

//depth of node
size_t StanTree::depth()
{
   if(!p) return 0; //no parents
   else return (1+p->depth());
}

//StanTree size
size_t StanTree::StanTreesize()
{
   if(l==0) return 1;  //if bottom node, StanTree size is 1
   else return (1+l->StanTreesize()+r->StanTreesize());
}

//node type
char StanTree::ntype()
{
   //t:top, b:bottom, n:no grandchildren, i:internal
   if(!p) return 't';
   if(!l) return 'b';
   if(!(l->l) && !(r->l)) return 'n';
   return 'i';
}

//print out StanTree(pc=true) or node(pc=false) information
void StanTree::pr(bool pc) 
{
   size_t d = depth();
   size_t id = nid();

   size_t pid;
   if(!p) pid=0; //parent of top node
   else pid = p->nid();

   std::string pad(2*d,' ');
   std::string sp(", ");
   if(pc && (ntype()=='t'))
      cout << "StanTree size: " << StanTreesize() << std::endl;
   cout << pad << "(id,parent): " << id << sp << pid;
   cout << sp << "(split_var,cut_val): " << split_var << sp << cut_val;
   cout << sp << "step_height: " << step_height;
   cout << sp << "type: " << ntype();
   cout << sp << "depth: " << depth();
   cout << sp << "pointer: " << this << std::endl;

   if(pc) {
      if(l) {
         l->pr(pc);
         r->pr(pc);
      }
   }
}

//kill children of  nog node nid
bool StanTree::death(size_t nid, double step_height)
{
   StanTree_p nb = getptr(nid);
   if(nb==0) {
      cout << "error in death, nid invalid\n";
      return false;
   }
   if(nb->isnog()) {
      delete nb->l;
      delete nb->r;
      nb->l=0;
      nb->r=0;
      nb->split_var=0;
      nb->cut_val=0;
      nb->step_height=step_height;
      return true;
   } else {
      cout << "error in death, node is not a nog node\n";
      return false;
   }
}

//is the node a nog node
bool StanTree::isnog() 
{
   bool isnog=true;
   if(l) {
      if(l->l || r->l) isnog=false; //one of the children has children.
   } else {
      isnog=false; //no children
   }
   return isnog;
}

size_t StanTree::nnogs() 
{
   if(!l) return 0; //bottom node
   if(l->l || r->l) { //not a nog
      return (l->nnogs() + r->nnogs());
   } else { //is a nog
      return 1;
   }
}

size_t StanTree::nbots() 
{
   if(l==0) { //if a bottom node
      return 1;
   } else {
      return l->nbots() + r->nbots();
   }
}

//get bottom nodes
void StanTree::getbots(npv& bv)
{
   if(l) { //have children
      l->getbots(bv);
      r->getbots(bv);
   } else {
      bv.push_back(this);
   }
}

//get nog nodes
void StanTree::getnogs(npv& nv)
{
   if(l) { //have children
      if((l->l) || (r->l)) {  //have grandchildren
         if(l->l) l->getnogs(nv);
         if(r->l) r->getnogs(nv);
      } else {
         nv.push_back(this);
      }
   }
}

//get all nodes
void StanTree::getnodes(npv& v)
{
   v.push_back(this);
   if(l) {
      l->getnodes(v);
      r->getnodes(v);
   }
}
void StanTree::getnodes(cnpv& v)  const
{
   v.push_back(this);
   if(l) {
      l->getnodes(v);
      r->getnodes(v);
   }
}

StanTree::StanTree_p StanTree::bn(double *x,xinfo& xi)
{
   if(l==0) return this; //no children
   if(x[split_var] < xi[split_var][cut_val]) {
      return l->bn(x,xi);
   } else {
      return r->bn(x,xi);
   }
}

//find region for a given variable
void StanTree::rg(size_t split_var, int* L, int* U)
{
   if(this->p==0)  {
      return;
   }
   if((this->p)->split_var == split_var) { //does my parent use split_var?
      if(this == p->l) { //am I left or right child
         if((int)(p->cut_val) <= (*U)) *U = (p->cut_val)-1;
         p->rg(split_var,L,U);
      } else {
         if((int)(p->cut_val) >= *L) *L = (p->cut_val)+1;
         p->rg(split_var,L,U);
      }
   } else {
      p->rg(split_var,L,U);
   }
}

//cut back to one node
void StanTree::tonull()
{
   size_t ts = StanTreesize();
   //loop invariant: ts>=1
   while(ts>1) { //if false ts=1
      npv nv;
      getnogs(nv);
      for(size_t i=0;i<nv.size();i++) {
         delete nv[i]->l;
         delete nv[i]->r;
         nv[i]->l=0;
         nv[i]->r=0;
      }
      ts = StanTreesize(); //make invariant true
   }
   step_height=0.0;
   split_var=0;cut_val=0;
   p=0;l=0;r=0;
}

//copy StanTree StanTree o to StanTree n
void StanTree::cp(StanTree_p n, StanTree_cp o)
//assume n has no children (so we don't have to kill them)
//recursion down
{
   if(n->l) {
      cout << "cp:error node has children\n";
      return;
   }

   n->step_height = o->step_height;
   n->split_var = o->split_var;
   n->cut_val = o->cut_val;

   if(o->l) { //if o has children
      n->l = new StanTree;
      (n->l)->p = n;
      cp(n->l,o->l);
      n->r = new StanTree;
      (n->r)->p = n;
      cp(n->r,o->r);
   }
}

//operators
StanTree& StanTree::operator=(const StanTree& rhs)
{
   if(&rhs != this) {
      tonull(); //kill left hand side (this)
      cp(this,&rhs); //copy right hand side to left hand side
   }
   return *this;
}

//functions
std::ostream& operator<<(std::ostream& os, const StanTree& t)
{
   StanTree::cnpv nds;
   t.getnodes(nds);
   os << nds.size() << std::endl;
   for(size_t i=0;i<nds.size();i++) {
      os << nds[i]->nid() << " ";
      os << nds[i]->GetSplitVar() << " ";
      os << nds[i]->GetCutVal() << " ";
      os << nds[i]->gettheta() << std::endl;
   }
   return os;
}
std::istream& operator>>(std::istream& is, StanTree& t)
{
   size_t tid,pid; //tid: id of current node, pid: parent's id
   std::map<size_t,StanTree::StanTree_p> pts;  //pointers to nodes indexed by node id
   size_t nn; //number of nodes

   t.tonull(); // obliterate old StanTree (if there)

   //read number of nodes----------
   is >> nn;
   if(!is) {
      //cout << ">> error: unable to read number of nodes" << endl;
      return is;
   }

   //read in vector of node information----------
   std::vector<node_info> nv(nn);
   for(size_t i=0;i!=nn;i++) {
      is >> nv[i].id >> nv[i].split_var >> nv[i].cut_val >> nv[i].step_height;
      if(!is) {
         //cout << ">> error: unable to read node info, on node  " << i+1 << endl;
         return is;
      }
   }
   //first node has to be the top one
   pts[1] = &t; //careful! this is not the first pts, it is pointer of id 1.
   t.setv(nv[0].split_var); t.setc(nv[0].cut_val); t.settheta(nv[0].step_height);
   t.p=0;

   //now loop through the rest of the nodes knowing parent is already there.
   for(size_t i=1;i!=nv.size();i++) {
      StanTree::StanTree_p np = new StanTree;
      np->split_var = nv[i].split_var; np->cut_val=nv[i].cut_val; np->step_height=nv[i].step_height;
      tid = nv[i].id;
      pts[tid] = np;
      pid = tid/2;
      // set pointers
      if(tid % 2 == 0) { //left child has even id
         pts[pid]->l = np;
      } else {
         pts[pid]->r = np;
      }
      np->p = pts[pid];
   }
   return is;
}

//add children to bot node *np
void StanTree::birthp(StanTree_p np,size_t split_var, size_t cut_val, double step_heightl, double step_heightr)
{
   StanTree_p l = new StanTree;
   l->step_height=step_heightl;
   StanTree_p r = new StanTree;
   r->step_height=step_heightr;
   np->l=l;
   np->r=r;
   np->split_var = split_var; np->cut_val=cut_val;
   l->p = np;
   r->p = np;
}

//kill children of  nog node *nb
void StanTree::deathp(StanTree_p nb, double step_height)
{
   delete nb->l;
   delete nb->r;
   nb->l=0;
   nb->r=0;
   nb->split_var=0;
   nb->cut_val=0;
   nb->step_height=step_height;
}

size_t StanTree::getbadcut(size_t split_var){
  StanTree_p par=this->getp();
  if(par->GetSplitVar()==split_var)
    return par->GetCutVal();
  else
    return par->getbadcut(split_var);
}

#ifndef NoRcpp   
// instead of returning y.test, let's return StanTrees
// this conveniently avoids the need for x.test
// loosely based on pr() 
// create an efficient list from a single StanTree
// StanTree2list calls itself recursively
Rcpp::List StanTree::StanTree2list(xinfo& xi, double center, double scale) {
  Rcpp::List res;

  // five possible scenarios
  if(l) { // StanTree has branches
    //double cut=xi[split_var][c];
    size_t var=split_var, cut=cut_val;

    var++; cut++; // increment from 0-based (C) to 1-based (R) array index

    if(l->l && r->l)         // two sub-StanTrees
      res=Rcpp::List::create(Rcpp::Named("var")=(int)var,
			     //Rcpp::Named("cut")=cut,
			     Rcpp::Named("cut")=(int)cut,
			     Rcpp::Named("type")=1,
			     Rcpp::Named("left")= l->StanTree2list(xi, center, scale),
			     Rcpp::Named("right")=r->StanTree2list(xi, center, scale));   
    else if(l->l && !(r->l)) // left sub-StanTree and right terminal
      res=Rcpp::List::create(Rcpp::Named("var")=(int)var,
			     //Rcpp::Named("cut")=cut,
			     Rcpp::Named("cut")=(int)cut,
			     Rcpp::Named("type")=2,
			     Rcpp::Named("left")= l->StanTree2list(xi, center, scale),
			     Rcpp::Named("right")=r->gettheta()*scale+center);    
    else if(!(l->l) && r->l) // left terminal and right sub-StanTree
      res=Rcpp::List::create(Rcpp::Named("var")=(int)var,
			     //Rcpp::Named("cut")=cut,
			     Rcpp::Named("cut")=(int)cut,
			     Rcpp::Named("type")=3,
			     Rcpp::Named("left")= l->gettheta()*scale+center,
			     Rcpp::Named("right")=r->StanTree2list(xi, center, scale));
    else                     // no sub-StanTrees 
      res=Rcpp::List::create(Rcpp::Named("var")=(int)var,
			     //Rcpp::Named("cut")=cut,
			     Rcpp::Named("cut")=(int)cut,
			     Rcpp::Named("type")=0,
			     Rcpp::Named("left")= l->gettheta()*scale+center,
			     Rcpp::Named("right")=r->gettheta()*scale+center);
  }
  else // no branches
    res=Rcpp::List::create(Rcpp::Named("var")=0, // var=0 means root
			   //Rcpp::Named("cut")=0.,
			   Rcpp::Named("cut")=0,
			   Rcpp::Named("type")=0,
			   Rcpp::Named("left") =step_height*scale+center,
			   Rcpp::Named("right")=step_height*scale+center);

  return res;
}

// for one StanTree, count the number of branches for each variable
Rcpp::IntegerVector StanTree::StanTree2count(size_t nvar) {
  Rcpp::IntegerVector res(nvar);

  if(l) { // StanTree branches
    res[split_var]++;
    
    if(l->l) res+=l->StanTree2count(nvar); // if left sub-StanTree
    if(r->l) res+=r->StanTree2count(nvar); // if right sub-StanTree
  } // else no branches and nothing to do

  return res;
}
#endif

