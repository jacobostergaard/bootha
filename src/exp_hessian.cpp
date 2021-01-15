// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <stdio.h>
#include <math.h>

using namespace std;
using namespace arma;
using namespace Rcpp;

// Declare functions

// [[Rcpp::export]]
arma::mat exp_hessian(arma::vec pp, double Tmax, arma::vec pars){
  
   /*
    *  pp:   observed point process
    *  Tmax: end time of the analyzed interval
    *  pars:  parameters: (mu, n, B)
    *          u:    background intensity
    *          y:     branching ratio
    *          b:     kernel parameter
    *  dt:   time resolution, if dt=0, then exact, but only for exponential kernel, dt>0 for generic kernels
    */
    
    mat out(pars.n_elem, pars.n_elem, fill::zeros);
    
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
      double Huu = -sum( pow(u+y*b*A,-2));
      double Huy = -sum( b*A % pow(u+y*b*A,-2) );
      double Hub = -sum( (y*A-y*b*B) % pow(u+y*b*A,-2) );
      
      double Hyy = -sum( pow(b*A,2) % pow(u+y*b*A,-2) );
      double Hyb = -sum( Tt%exp(-b*Tt) ) + sum( (A-b*B)/(u+y*b*A) ) - sum( y*b*A%(A-b*B)%pow(u+y*b*A,-2) );
      
      double Hbb = y*sum( exp(-b*Tt)%pow(Tt,2) )+ sum( (y*b*C-2*y*B)/(u+y*b*A) ) - sum( pow(y*A-y*b*B,2) % pow(u+y*b*A,-2) );
      
      out(0,0) = Huu;
      out(0,1) = Huy;
      out(0,2) = Hub;
      
      out(1,0) = Huy;
      out(1,1) = Hyy;
      out(1,2) = Hyb;
      
      out(2,0) = Hub;
      out(2,1) = Hyb;
      out(2,2) = Hbb;
    
      return out;
}
