// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <stdio.h>
#include <math.h>

using namespace std;
using namespace arma;
using namespace Rcpp;

// Declare functions
arma::vec exp_kernel(arma::vec x, double a, double b); 
arma::vec pow_kernel(arma::vec x, double a, double b, double k);
  

// [[Rcpp::export]]
arma::vec intensity(arma::vec x, arma::vec pp, arma::vec pars){
      
    /*
     *  x:     evaluation point
     *  pp:    observed point process
     *  pars:  parameters: (mu, n, B) or (mu, a, b, k) , ie exp or powerlaw
     *          mu:    background intensity
     *          n:     branching ratio
     *          B:     kernel parameter
     *          
     *          a, b, k are parameters for the powerlaw kernel
     *
     */

    int nx   = x.n_elem;       // number of evaluation points
    int np   = pars.n_elem;    // number of parameters (3 pars=> exp kernel, 4 pars => powerlaw)
    vec out  = zeros<vec>(nx);
    vec tmp;
    uvec idx;

    for(int i=0; i<nx; i++){
      idx = find( pp < x(i) );               // which are events prior to current evaluation point
      if(idx.n_elem>0){
        tmp = x(i)-pp(idx);                    // lagged times to use in kernel evaluation  
        if(np==3){ // use exponential kernel
          out(i) = pars(0)+sum( exp_kernel(tmp, pars(1)*pars(2), pars(2))  );  
        } else if(np==4){ // use powerlaw kernel
          out(i) = pars(0)+sum( pow_kernel(tmp, pars(1), pars(2), pars(3))  );
        }
      } else{ 
        out(i) = pars(0);
      }
    }

  return out;
}


// [[Rcpp::export]]
arma::vec integrated_intensity(arma::vec x, arma::vec pp, arma::vec pars){
    /*
     *  x:     evaluation point
     *  pp:    observed point process
     *  pars:  parameters: (mu, n, B)
     *          mu:    background intensity
     *          n:     branching ratio
     *          B:     kernel parameter
     *
     */
    
  int nx = x.n_elem;
  vec out = zeros<vec>(nx);
  vec tmp; 
  
  tmp = intensity(x, pp, pars);               // find intensity at x-points
  tmp = .5*tmp.tail(nx-1)+.5*tmp.head(nx-1);  // integrate using over/undersums
  
  vec dt  = x.tail(nx-1)-x.head(nx-1);        // set the time scale resolution
  
  out.tail(nx-1) = cumsum(tmp%dt);            // integration result
  
  return out;
}