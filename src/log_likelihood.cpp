// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <stdio.h>
#include <math.h>

using namespace std;
using namespace arma;
using namespace Rcpp;

// Declare functions
arma::vec intensity(arma::vec x, arma::vec pp, arma::vec pars);                              // conditional intensity
arma::vec integrated_intensity(arma::vec x, arma::vec pp, arma::vec pars);                              // integreated intensity

// [[Rcpp::export]]
double log_lik(arma::vec pp, double Tmax, arma::vec pars,  double dt = 0){

  /*
   *  pp:   observed point process
   *  Tmax: end time of the analyzed interval
   *  pars:  parameters: (mu, n, B)
   *          mu:    background intensity
   *          n:     branching ratio
   *          B:     kernel parameter
   *  dt:   time resolution, if dt=0, then exact, but only for exponential kernel, dt>0 for generic kernels
   */

  double out;
  if(dt == 0){
      // Fast likelihood evaluation with exp kernel:
      
        double mu = pars(0);
        double n = pars(1);
        double B = pars(2);
      
        double tmp1 = -mu*Tmax;
        double tmp2 = n*sum( exp( -B*(Tmax-pp) )-1 );
    
        int np   = pp.n_elem;
        vec A    = zeros<vec>(np);
    
        for(int i=1; i<np; i++){
          A(i) = exp( -B*(pp(i)-pp(i-1)) )*( 1+A(i-1) );
        }
        double tmp3 = sum(log(mu+ n*B*A));
    
        out  = tmp1+tmp2+tmp3;
  } else {
    // Generic likelihood evalutation (slow due to integration step):
    // double dt = 0.1;
    vec t = regspace(0,dt,Tmax);
    double tmp1 = -sum(intensity(t,pp,pars))*dt;
    double tmp2 = sum(log(intensity(pp, pp, pars)));
    out = tmp1+tmp2;
  }
    
    return out;
}
