#include <Rcpp.h>
#include <iostream>
#include <string>
#include <random>
#include <boost/math/distributions/beta.hpp>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector timesTwo(NumericVector x) {
  return x * 2;
}



