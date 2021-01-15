// #include <Rcpp.h>
// using namespace Rcpp;
#include <RcppArmadillo.h>
#include <stdio.h>
#include <math.h>

using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]
using namespace std;
using namespace arma;

arma::vec exp_kernel(arma::vec x, double a, double b){
  
  vec out = a*exp(-b*(x));
    
  return out;
}
  
arma::vec pow_kernel(arma::vec x, double a, double b, double k){
    
  vec out = (a*b)/pow(1+b*x,1+k);
    
  return out;
}
