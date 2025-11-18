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
  virtual size_t discrete() = 0;  // Discrete (categorical) distribution
  virtual void SetInclusionWeights(std::vector<double>& _wts) = 0;
  virtual ~Random() {}
};

// Abstract random number generator based on C++ <random>
class arn : public Random {
  typedef std::default_random_engine genD;
  typedef std::normal_distribution<double> norD;
  typedef std::uniform_real_distribution<double> uniD;
  typedef std::chi_squared_distribution<double> chiD;
  typedef std::gamma_distribution<double> gamD;
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
  
  virtual size_t discrete() {
    disD dis(wts.begin(), wts.end());
    return dis(gen);
  }

  virtual void SetInclusionWeights(std::vector<double>& _wts) {
    double smw = 0.;
    wts.clear();
    for (size_t j = 0; j < _wts.size(); j++) smw += _wts[j];
    for (size_t j = 0; j < _wts.size(); j++) wts.push_back(_wts[j] / smw);
  }

private:
  std::vector<double> wts;
  genD gen;
  norD nor;
  uniD uni;
  chiD chi;
  gamD gam;
};

#endif // GUARD_Random_h