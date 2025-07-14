#ifndef GUARD_Random_h
#define GUARD_Random_h

#include <vector>
#include <random>
#include <cmath>

// Pure virtual base class for random numbers
class Random {
public:
  Random() {}
  virtual double normal() = 0;    // Standard normal
  virtual double uniform() = 0;   // Uniform(0,1)
  virtual double chi_square(double df) = 0; // Chi-square
  virtual double exp() = 0;       // Exponential
  virtual double log_gamma(double shape) = 0;
  virtual double gamma(double shape, double rate) = 0;
  virtual double inv_gamma(double shape, double rate) = 0;
  virtual double beta(double a, double b) = 0;
  virtual size_t discrete() = 0;  // Discrete (categorical) distribution
  virtual size_t geometric(double p) = 0;   // Geometric distribution
  virtual double inverse_gaussian(double mu, double lambda) = 0; // Inverse Gaussian
  virtual void SetInclusionWeights(std::vector<double>& _wts) = 0;
  virtual std::vector<double> log_dirichlet(std::vector<double>& alpha) = 0;
  virtual ~Random() {}
};

// Abstract random number generator based on C++ <random>
class arn : public Random {
  typedef std::default_random_engine genD;
  typedef std::normal_distribution<double> norD;
  typedef std::uniform_real_distribution<double> uniD;
  typedef std::chi_squared_distribution<double> chiD;
  typedef std::gamma_distribution<double> gamD;
  typedef std::geometric_distribution<int> geoD;
  typedef std::discrete_distribution<int> disD;
  
public:
  // Constructor
  arn() : gen(std::random_device{}()) {}
  arn(unsigned int n1, unsigned int n2) : gen(n1 + n2) {}
  
  // Virtual destructor
  virtual ~arn() {}
  
  // Implementations of pure virtual functions
  virtual double normal() { return nor(gen); }
  
  virtual double uniform() { return uni(gen); }
  
  virtual double chi_square(double df) {
    chi = chiD(df);
    return chi(gen);
  }
  
  virtual double exp() { return -std::log(this->uniform()); }
  
  virtual double log_gamma(double shape) {
    gam = gamD(shape + 1., 1.);
    double y = std::log(gam(gen)), z = std::log(this->uniform()) / shape;
    return y + z;
  }
  
  virtual double gamma(double shape, double rate) {
    if (shape < 0.01) {
      return std::exp(this->log_gamma(shape)) / rate;
    } else {
      gam = gamD(shape, 1.);
      return (gam(gen)) / rate;
    }
  }
  
  virtual double inv_gamma(double a, double b) {
    double gamma_sample = this->gamma(a, b);
    return 1.0 / gamma_sample;
  }
  
  virtual double beta(double a, double b) {
    double x1 = this->gamma(a, 1.), x2 = this->gamma(b, 1.);
    return x1 / (x1 + x2);
  }
  
  virtual size_t discrete() {
    disD dis(wts.begin(), wts.end());
    return dis(gen);
  }
  
  virtual size_t geometric(double p) {
    geo = geoD(p);
    return geo(gen);
  }
  
  virtual double inverse_gaussian(double mu, double lambda) {
    std::normal_distribution<> normal_dist(0.0, 1.0);
    std::uniform_real_distribution<> uniform_dist(0.0, 1.0);
    
    double Z = normal_dist(gen);  // Step 1: Standard normal random variable
    double Y = Z * Z;             // Step 2: Y = Z^2
    
    double mu_squared = mu * mu;
    double X = mu + (mu_squared * Y) / (2 * lambda) 
      - (mu / (2 * lambda)) * std::sqrt(4 * mu * lambda * Y + mu_squared * Y * Y);
    
    // Step 4: Acceptance-rejection criterion
    double U = uniform_dist(gen);
    if (U <= mu / (mu + X)) {
      return X;
    } else {
      return mu_squared / X;
    }
  }
  
  virtual void SetInclusionWeights(std::vector<double>& _wts) {
    double smw = 0.;
    wts.clear();
    for (size_t j = 0; j < _wts.size(); j++) smw += _wts[j];
    for (size_t j = 0; j < _wts.size(); j++) wts.push_back(_wts[j] / smw);
  }
  
  virtual std::vector<double> log_dirichlet(std::vector<double>& alpha) {
    size_t k = alpha.size();
    std::vector<double> draw(k);
    double lse;
    for (size_t j = 0; j < k; j++) draw[j] = this->log_gamma(alpha[j]);
    
    double max_elem = *std::max_element(draw.begin(), draw.end());
    double sum = 0.0;
    for (double x : draw) {
      sum += std::exp(x - max_elem);
    }
    
    lse = max_elem + std::log(sum);
    for (size_t j = 0; j < k; j++) {
      draw[j] -= lse;
    }
    return draw;
  }
  
private:
  std::vector<double> wts;
  genD gen;
  norD nor;
  uniD uni;
  chiD chi;
  gamD gam;
  geoD geo;
};

#endif
