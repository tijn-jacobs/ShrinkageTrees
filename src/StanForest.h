
#ifndef GUARD_StanForest_h
#define GUARD_StanForest_h

#include <ctime>

#include "StanTree.h"
#include "StanTreeFunctions.h"
#include "Info.h"
#include "StanForestFunctions.h"
#include "bd.h"

class StanForest {
public:
   
   //friends
   friend bool bd(StanTree& x, xinfo& xi, dinfo& di, pinfo& pi, double sigma,
		  std::vector<size_t>& nv, std::vector<double>& pv, bool aug, Random& random);
   
   //constructor/destructor
   // StanForest();
   StanForest(size_t m);
   // StanForest(const StanForest&);
   ~StanForest();
   
   //operators
   StanForest& operator=(const StanForest&);
   
   //get,set
   size_t getm() {return m;}
   void setm(size_t m);
   void setdata(size_t p, size_t n, double *x, double *y, size_t nc=100);
   void setdata(size_t p, size_t n, double *x, double *y, int* nc);
   void setpi(pinfo& pi) {this->pi = pi;}
   void setprior(double base, double power, double eta)
      {pi.base=base; pi.power = power; pi.eta=eta;}
   void setdart(double _a, double _b, double _rho, bool _aug, bool _dart, 
		double _theta=0.) {
     this->a=_a; this->b=_b; this->rho=_rho; this->aug=_aug; 
     this->dart=_dart;  
     if(_theta==0.){
       this->const_theta=false;
       this->theta=1.;
     }
     else {
       this->const_theta=true;
       this->theta=_theta;
     }
}
   void startdart() {this->dartOn=!(this->dartOn);}
   StanTree& getStanTree(size_t i ) { return t[i];}
   xinfo& getxinfo() {return xi;}
   void setxinfo(xinfo& _xi);
   std::vector<size_t>& getnv() {return nv;}
   std::vector<double>& getpv() {return pv;}
   double gettheta() {return theta;}
   
   //public methods
   void birth(size_t i, size_t nid,size_t v, size_t c, double ml, double mr)
         {t[i].birth(nid,v,c,ml,mr);}
   void death(size_t i,size_t nid, double mu)
         {t[i].death(nid,mu);}
   void pr();
   void tonull() {for(size_t i=0;i!=t.size();i++) t[i].tonull();}
   void predict(size_t p, size_t n, double *x, double *fp);
   bool draw(double sigma, Random& random, bool* accept);
   double f(size_t i) {return allfit[i];}
   double* GetAllFit() { return allfit; }
  void UpdateGlobalScaleParameters(string prior_type,
                                            double global_parameter,
                                            double& storage_eta, // Store the updated eta at this location
                                            Random& random);
  void UpdateHalfCauchyScale(double global_parameter, double& storage_eta, Random& random);

protected:
   size_t m;  //number of StanTrees
   std::vector<StanTree> t; //the StanTrees
   pinfo pi; //prior and mcmc info
   //data
   size_t p,n; //x has dim p, n obserations
   double *x,*y;  //x is column stack, pxn
   xinfo xi; //cutpoint info
   //working
   double *allfit; //if the data is set, should be f(x)
   double *r;
   double *ftemp;
   dinfo di;
   bool dart,dartOn,aug,const_theta;
   double a,b,rho,theta;
   std::vector<size_t> nv;
   std::vector<double> pv, lpv;
};

#endif
