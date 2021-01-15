// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <stdio.h>
#include <math.h>

using namespace std;
using namespace arma;
using namespace Rcpp;

// Declare functions

// [[Rcpp::export]]
arma::vec exp_gradient(arma::vec pp, double Tmax, arma::vec pars){
  
  /*
   *  pp:   observed point process
   *  Tmax: end time of the analyzed interval
   *  pars:  parameters: (mu, n, B)
   *          u:    background intensity
   *          y:     branching ratio
   *          b:     kernel parameter
   *  dt:   time resolution, if dt=0, then exact, but only for exponential kernel, dt>0 for generic kernels
   */
    
    vec out = zeros<vec>(pars.n_elem);
    
    // Fast likelihood evaluation with exp kernel:
  
    double u = pars(0);
    double y = pars(1);
    double b = pars(2);
    
    int n  = pp.n_elem;
    vec A  = zeros<vec>(n);
    vec B  = zeros<vec>(n);
    vec C  = zeros<vec>(n);
    
    double ti;
    
    for(int i=1; i<n; i++){
      ti   = pp(i)-pp(i-1);
    
      A(i) = exp( -b*ti )*( 1+A(i-1) );
      
      B(i) = ti*A(i)+exp(-b*ti)*B(i-1);
      C(i) = ti*(B(i)+exp(-b*ti)*B(i-1))+exp(-b*ti)*C(i-1);
      
    }
    
    vec Tt = Tmax-pp;
    double dldu = -Tmax+sum( 1/(u+y*b*A) );
    double dldy = sum(  exp(-b*Tt)-1 + b*A/(u+y*b*A)  );
    double dldb = -y*sum( exp(-b*Tt)%Tt  ) + sum( (y*A-y*b*B)/(u+y*b*A) );


    out(0) = dldu;
    out(1) = dldy;
    out(2) = dldb;
      
  return out;
}
