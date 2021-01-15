// #include <Rcpp.h>
#include <RcppArmadillo.h>
#include <stdio.h>
#include <math.h>

using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]
using namespace std;
using namespace arma;

// Declare functions
arma::vec interpol_aux(arma::mat f, arma::vec x, bool flat=false);
arma::vec integrated_intensity(arma::vec x, arma::vec pp, arma::vec pars);



double g(double x, arma::vec pars, double tk, double Sk){
  // Used for fast evaluation of rib for exponential kernel
  double out = pars(0)*(x-tk)+pars(1)*(1-exp(-pars(2)*(x-tk)))*Sk;
  
  return out;
}

double dg(double x, arma::vec pars, double tk, double Sk){
  // Derivative of f, used for fast evaluation of rib for exponential kernel (to find 0's with Newtons algorithm)
  double out = pars(0)+pars(1)*Sk*exp(-pars(2)*(x-tk));
  return out;
}

// [[Rcpp::export]]
arma::vec fib(arma::mat L, arma::vec v, bool flat=false) {
// arma::vec fib(arma::vec pars, arma::vec pp, arma::vec s, double dt = 1, bool flat=false) {

  /*
   *    Inputs:
   *      pars  parameters for the intensity function
   *      pp    observed point process, used to find the observed integrated intensity, given parameters
   *      s     set of waiting times, either unit exponential Exp(1) or from data
   *      dt    time resolution for integration
   *      
   */
  vec out;
  vec s = cumsum(v);
  mat h = L;
   
  h.col(0) = L.col(1);
  h.col(1) = L.col(0);
  // vec x   = regspace(0,dt,1.5*max(pp)); 
  // vec tmp = integrated_intensity(x, pp, pars);
  // mat g   = join_rows(tmp,x);
  
  out = interpol_aux(h,s,flat);
  
  return out;
  
}



// [[Rcpp::export]]
arma::vec exp_fib(arma::mat L, arma::vec v, arma::vec pars, arma::vec pp, bool flat=false) {
  
  
  // vec x   = regspace(0,dt,1.5*max(pp)); 
  // vec L = integrated_intensity(x, pp, pars);
  // mat h   = join_rows(x,L);
  

  vec tmp = interpol_aux(L,pp,flat); // time transformed observed events
  
  int n = v.n_elem;
  
  vec out = zeros<vec>(n);
  vec s = cumsum(v);
  
  double tk;
  double Sk;
  double sk;
  
  vec idx;
  int k;
  
  double x0;
  double x1; 
  double xd; 
  double tol = 1e-6;
  int it;
  int niter = 1e3;
  
  for(int i=0; i<n; i++){
    
    // idx = find( s(i) > tmp );
    if( any(s(i) > tmp) ){
      k  = max( find( s(i) > tmp ) );
      tk = pp(k);
      Sk = sum( exp( -pars(2)*(tk-pp.subvec(0,k) ) ));
      sk = tmp(k);
    } else{
      tk = 0;
      Sk = 0;
      sk = 0;
    }
  
    if(i>0){
      x0 = out(i-1);
    } else{
      x0 = 0;
    }
    
    // solve f(x)=0 by Newtons algorithm:
    xd = 1;
    it = 0;

    while(xd > tol & it < niter){
      
      x1 = x0-(g(x0, pars, tk, Sk)+sk-s(i))/dg(x0,pars,tk,Sk);
      xd = x1-x0;
      x0 = x1;
      it = it+1;
    }
        // update observed times:
      out(i) = x1;
  } 
    
  return out;
}



// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R

library(bootha)
library(misc)
misc::clean_up()
set.seed(1234)

plotlib = "~/Library/Mobile Documents/com~apple~CloudDocs/GitHub/Source/R/packages/bootha/plots/"
datalib = "~/Library/Mobile Documents/com~apple~CloudDocs/GitHub/Source/R/packages/bootha/data/"
load(paste0(datalib,"boot_embrechts5.Rda"))
load(paste0(datalib,"dowjones.Rda"))

pp   = dowjones$pp.neg
pars = dowjones$pars
Tmax = dowjones$Tmax

mle  = mle(pp, Tmax, pars, dt=0, grad = function(x) exp_gradient(pp, Tmax, x), hessian=FALSE)
n  = 50
v  = rexp(n)

s  = cumsum(v)

x = seq(0,Tmax, .01)
L = cbind(x,integrated_intensity(x, pp, mle))

system.time({
  num_fib = fib(L, s)  
})

system.time({
  exp_fib = exp_fib(L, s, mle, pp)
})

sum((exp_fib-num_fib)^2)
plot(exp_fib, num_fib);abline(0,1)

  */
