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

// What does this function do?
// It finds the bottom node for a given x. It uses the split rules in the tree 
// to find the bottom node. It returns a pointer to the bottom node associated with x.
StanTree::StanTree_p StanTree::bn(double *x, xinfo& xi) {
   if(l==0) return this; //no children
   if(std::isnan(x[split_var])) {
    // This is for missing data!
    } else {
      if(x[split_var] < xi[split_var][cut_val]) {
        return l->bn(x,xi);
      } else {
        return r->bn(x,xi);
      }
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
