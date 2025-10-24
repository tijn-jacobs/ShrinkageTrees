// Timing.h
#pragma once
#include <chrono>
#include <unordered_map>
#include <string>

struct ProfAgg {
  static std::unordered_map<std::string, double>& bucket() {
    static std::unordered_map<std::string, double> b;
    return b;
  }
};

struct ScopedTimer {
  using clock = std::chrono::high_resolution_clock;
  const std::string name;
  std::chrono::time_point<clock> t0;
  ScopedTimer(const std::string& n) : name(n), t0(clock::now()) {}
  ~ScopedTimer() {
    auto t1 = clock::now();
    double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    ProfAgg::bucket()[name] += ms;
  }
};

// [[Rcpp::export]]
Rcpp::NumericVector get_cpp_profile() {
  Rcpp::NumericVector out;
  for (auto& kv : ProfAgg::bucket()) out[kv.first] = kv.second;
  return out;
}

// [[Rcpp::export]]
void reset_cpp_profile() { ProfAgg::bucket().clear(); }
