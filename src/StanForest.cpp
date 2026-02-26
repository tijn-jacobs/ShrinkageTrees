#include "StanForest.h"


// Constructor
StanForest::StanForest(size_t im):m(im),t(m),pi(),p(0),n(0),x(0),y(0),xi(),allfit(0),r(0),ftemp(0),di(),dartOn(false) {}

// Destructor
StanForest::~StanForest()
{
   if(allfit) delete[] allfit;
   if(r) delete[] r;
   if(ftemp) delete[] ftemp;
}

// Operators
StanForest& StanForest::operator=(const StanForest& rhs)
{
   if(&rhs != this) {

      this->t = rhs.t;
      this->m = t.size();

      this->pi = rhs.pi;

      p=0;n=0;x=0;y=0;
      xi.clear();

      if(allfit) {delete[] allfit; allfit=0;}
      if(r) {delete[] r; r=0;}
      if(ftemp) {delete[] ftemp; ftemp=0;}

   }
   return *this;
}

//get,set
void StanForest::setm(size_t m)
{
   t.resize(m);
   this->m = t.size();

   if(allfit && (xi.size()==p)) predict(p,n,x,allfit);
}


void StanForest::setxinfo(xinfo& _xi)
{
   size_t p=_xi.size();
   xi.resize(p);
   for(size_t i=0;i<p;i++) {
     size_t nc=_xi[i].size();
      xi[i].resize(nc);
      for(size_t j=0;j<nc;j++) xi[i][j] = _xi[i][j];
   }
}

void StanForest::setdata(size_t p, size_t n, double *x, double *y, size_t numcut)
{
  int* nc = new int[p];
  for(size_t i=0; i<p; ++i) nc[i]=numcut;
  this->setdata(p, n, x, y, nc);
  delete [] nc;
}

void StanForest::setdata(size_t p, size_t n, double *x, double *y, int *nc)
{
   this->p=p; this->n=n; this->x=x; this->y=y;
   if(xi.size()==0) makexinfo(p,n,&x[0],xi,nc); // This should be always the case; CHECK!

   if(allfit) delete[] allfit;
   allfit = new double[n];
   predict(p,n,x,allfit);

   if(r) delete[] r;
   r = new double[n];

   if(ftemp) delete[] ftemp;
   ftemp = new double[n];

   di.n=n; di.p=p; di.x = &x[0]; di.y=r;

   nv.clear();
   pv.clear();
   for(size_t j=0;j<p;j++){
     nv.push_back(0);
     pv.push_back(1/(double)p);
   }
}

void StanForest::predict(size_t p, size_t n, double *x, double *fp)
//uses: m,t,xi
{
   double *fptemp = new double[n];

   for(size_t j=0;j<n;j++) fp[j]=0.0;
   for(size_t j=0;j<m;j++) {
      fit(t[j],xi,p,n,x,fptemp);
      for(size_t k=0;k<n;k++) fp[k] += fptemp[k];
   }

   delete[] fptemp;
}

bool StanForest::draw(double sigma, Random& random, bool* accept) {
   for(size_t j = 0; j < m; j++) {

      // remove old contribution
      fit(t[j], xi, p, n, x, ftemp);
      for(size_t k = 0; k < n; k++) {
         allfit[k] = allfit[k] - ftemp[k];
         r[k] = y[k] - allfit[k];
      }

      // store acceptance for this tree
      accept[j] = bd(t[j], xi, di, pi, sigma, nv, pv, aug, random);

      // draw new mu
      drmu(t[j], xi, di, pi, sigma, random);

      // add updated contribution
      fit(t[j], xi, p, n, x, ftemp);
      for(size_t k = 0; k < n; k++) {
         allfit[k] += ftemp[k];
      }
   }
   
   if(dartOn) {
     draw_s(nv, lpv, theta, random);
     draw_theta0(const_theta, theta, lpv, a, b, rho, random);
     for(size_t j = 0; j < p; j++)
         pv[j] = ::exp(lpv[j]);
   }

   return true;  // or return nothing if function is void
}


//public functions
void StanForest::pr() //print to screen
{
   cout << "*****StanForest object:\n";
   cout << "m: " << m << std::endl;
   cout << "t[0]:\n " << t[0] << std::endl;
   cout << "t[m-1]:\n " << t[m-1] << std::endl;
   cout << "prior and mcmc info:\n";
   if(dart){
     cout << "*****dart prior (On):\n";
     cout << "a: " << a << std::endl;
     cout << "b: " << b << std::endl;
     cout << "rho: " << rho << std::endl;
     cout << "augmentation: " << aug << std::endl;
   }
   else cout << "*****dart prior (Off):\n";
   if(p) cout << "data set: n,p: " << n << ", " << p << std::endl;
   else cout << "data not set\n";
}

void StanForest::UpdateGlobalScaleParameters(string prior_type,
                                            double global_parameter,
                                            double& storage_eta, // Store the updated eta at this location
                                            Random& random) {

  if (prior_type == "standard-halfcauchy") {
    this->UpdateHalfCauchyScale(global_parameter, storage_eta, random);
  } else if (prior_type == "standard-halfnormal") {
    return;
  } else { // Implement forest-wide horseshoe shrinkage update here
    return;  
  }  
}


void StanForest::UpdateHalfCauchyScale(double global_parameter, double& storage_eta, Random& random) {
    
  double sum_sq = 0.0;
  double leaf_count = 0.0;

  // Loop over all trees
  for (size_t j = 0; j < m; j++) {

      StanTree::npv leaves;
      t[j].getbots(leaves); // bottom nodes

      for (auto* leaf : leaves) {

          double h = leaf->gettheta();  // leaf parameter (scalar)

          sum_sq += (h * h) / (pi.eta * pi.eta);
          leaf_count += 1.0;
      }
  }

  double shape = 0.5 * (1.0 + leaf_count);
  double rate  = 0.5 * (1.0 + sum_sq);

  // Horseshoe auxiliary update
  double aux = random.gamma(shape, 1.0) / rate;

  // Update the global scale parameter
  pi.eta = global_parameter / (std::sqrt(aux) * std::sqrt((double) m));
  storage_eta = pi.eta;
}




