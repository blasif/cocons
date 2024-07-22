#define BOOST_DISABLE_ASSERTS
#ifndef AUXFUNCTIONS_H
#define AUXFUNCTIONS_H

#include <cstdlib>
#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <limits>
#include <boost/math/special_functions/bessel.hpp>

inline double Pexpfma_new(const Rcpp::NumericVector& b, const Rcpp::NumericVector& x){
  double tmp_sum = 0.0;
  for (int i = 0; i < x.length(); ++i) {
    tmp_sum = std::fma(x(i), b(i), tmp_sum);
  }
  return (1 / std::exp(- 1 * tmp_sum));
}

inline double Pexpfma_new_smoothness(const Rcpp::NumericVector& b, 
                                     const Rcpp::NumericVector& x, 
                                     double min_v, double max_v){
  double tmp_sum = 0.0;
  for (int i = 0; i < x.length(); ++i) {
    tmp_sum = std::fma(x(i), b(i), tmp_sum);
  }
  return ((max_v - min_v) / (1 + std::exp(- 1 * tmp_sum)) + min_v);
}

inline double Pexpfma_new_smoothness_two(const Rcpp::NumericVector& b, 
                                     const Rcpp::NumericVector& x, 
                                     double min_v, double max_v){
  double tmp_sum = 0.0;
  for (int i = 1; i < x.length(); ++i) {
    tmp_sum = std::fma(x(i), b(i), tmp_sum);
  }
  return ((max_v - min_v) / (1 + std::exp(- 1 * tmp_sum)) + min_v);
}

inline double newinvlogitfma(const Rcpp::NumericVector& b, 
                             const Rcpp::NumericVector& x){
  double tmp_sum = 0.0;
  for (int i = 0; i < x.length(); ++i) {
    tmp_sum = std::fma(x(i), b(i), tmp_sum);
  }
  return (M_PI / (1 + (std::exp( - 1 * tmp_sum)))) ;
}

inline double kahan(double a, double b, double c, double d) {
  double cd = c * d;
  double error = std::fma(c, d, -cd);
  double result = std::fma(a, b, -cd);
  return (result - error);
}

#endif // AUX_FUNCTIONS_H
