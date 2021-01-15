// #include <Rcpp.h>
#include <RcppArmadillo.h>
#include <stdio.h>
#include <math.h>
  
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
using namespace std;
using namespace arma;

// Declare functions
arma::vec exp_kernel(arma::vec x, double a, double b); 
arma::vec pow_kernel(arma::vec x, double a, double b, double k);
arma::vec interpol_aux(arma::mat f, arma::vec x, bool flat=false);
arma::vec integrated_intensity(arma::vec x, arma::vec pp, arma::vec pars);
  

double f(double x, arma::vec pars, double Sm, double v){
  // Used for fast evaluation of rib for exponential kernel
  double tmp = exp(-pars(2)*x);
  double out = pars(0)*x+pars(1)*(1-tmp)*Sm-v;
  
  return out;
}

double df(double x, arma::vec pars, double Sm){
  // Derivative of f, used for fast evaluation of rib for exponential kernel (to find 0's with Newtons algorithm)
  double out = pars(0)+pars(1)*pars(2)*exp(-pars(2)*x)*Sm;
  return out;
}
  
// [[Rcpp::export]]
arma::vec rib(arma::vec pars, arma::vec v, double dt = 1, bool flat=false) {
  
  /*
   *    Inputs:
   *      pars  parameters for the use kernel representation
   *      s     Set of waiting times, either unit exponential Exp(1) or from data
   *      dt    time resolution for integration
   *      flat  interpolation parameter, should be set to false
   */
    double mu_tmp = 1/pars(0);
    
    vec out = zeros<vec>(v.n_elem);
    vec s = cumsum(v);
    vec tmp;
    vec t;
    mat f;
    double T_tmp;
    double int_max;
    out(0) = s(0)/pars(0); // first event is simple due to linear integrated intensity (mu only)
    
    for(int i=1; i<s.n_elem; i++){
      int k=1;
      int_max = -1;
      T_tmp   = s(i)*mu_tmp+out(i-1); // estimated Tmax, based on parameters,  
      
      while(int_max <= 0){ // continue untill 0 can be interpolated
        t       = regspace(out(i-1),dt,k*T_tmp); // times k to try and reach sufficient Tmax!
        tmp     = integrated_intensity(t,out.subvec(0,i-1),pars)+s(i-1)-s(i);
        int_max = max(tmp);
        k       = k+1; // increase time interval to try and capture 0 in the S-time
      }

      f     = join_rows(tmp,t);
      out(i) = as_scalar( interpol_aux(f, zeros<vec>(1), flat) );
    }
    
    
  return out;
}


// [[Rcpp::export]]
arma::vec exp_rib(arma::vec pars, arma::vec v) {
  // fast rib evaluation using exact formulas and solving f(x)=0 by Newtons algorithm
  
  // initialize 
  int n  = v.n_elem;
  vec w  = zeros<vec>(n);
  w(0)   = v(0)/pars(0); // first waiting time depends only on baseline
  // vec v  = s;
  // v.subvec(1,n) = diff(s);
  
  double Sm = 1;
  double x0;
  double x1; 
  double xd; 
  double tol = 1e-6;
  int it;
  int niter = 100;
  
  for(int i = 1; i<n; i++){
    
    // solve f(x)=0 by Newtons algorithm:
    x0 = 0;
    xd = 1;
    it = 0;
    while(xd > tol & it < niter){
      x1 = x0-f(x0,pars,Sm,v(i))/df(x0,pars,Sm);
      xd = x1-x0;
      x0 = x1;
      it = it+1;
    }
    // update waiting time w and Sm
    w(i) = x0;
    Sm = exp( -pars(2)*w(i) )*Sm+1;
  }
  
  vec out = cumsum(w);
//   for(i in 2:n){
//    
//     while(xdiff > tol & iter < 10){
//       x1 = x0-f(x0,O,Sm,v[i])/dfdx(x0,O,Sm)
//       xdiff = x1-x0
//       x0 = x1
//       iter = iter+1
// # cat("\niter",iter,"x0:",x0)
//     }
// # w[i] = uniroot(f, c(0,100))$root
//     w[i] = x0
//       Sm = exp(-O[3]*w[i])*Sm+1
//   }
//   
//   out = cumsum(w)
  
  
  return out;
}






// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//
  
  /*** R
  misc::clean_up()
  
  v  = rexp(5) # sample events in S-time
  pars = c(0.018, 0.74, 0.021)


  # out = s[1]/pars[1]
  # mu.tmp = (1-pars[2])/pars[1]
  # for(i in 2:length(s)){
  #   T.tmp = s[i]*mu.tmp+out[i-1]
  #   k=1
  #   Z.int = -1
  #   while(max(Z.int)<=0){
  #     print(k)
  #     x.tmp = seq(out[i-1],k*T.tmp,.5)
  #     Z.tmp = integrated_intensity(x.tmp,out[1:(i-1)], pars)+s[i-1]
  #     Z.int = Z.tmp-s[i]
  #     k=k+1
  #   }
  #   out[i] = interpol(cbind(x.tmp,Z.int), y=0)
  # }
  # 
  # cbind( out,rib(pars,s, dt=.5) )
  
  
  # x = seq(0,max(out),.1)
  # Z = integrated_intensity(x, out, pars)
  # plot(x, Z, type='l', ylim=c(0,20))
  # abline(v=out, h=s, lty=3)
  # for(i in 2:length(s)){
  #   T.tmp = s[i]*mu.tmp+out[i-1]
  #   x.tmp = seq(out[i-1],T.tmp,.5)
  #   Z.tmp = integrated_intensity(x.tmp,out[1:(i-1)], pars)+s[i-1]
  #   points(x.tmp, Z.tmp, col='red', pch=16,cex=.25)
  #   Sys.sleep(1)
  # }

  # x1 = seq(0,3*T.tmp,1)
  # Z1 = integrated_intensity(x1,out[1:(i-1)], pars)-s[i]
  
  # plot(x1,Z1, type='l', ylim=c(-2,2))
  # abline(h=0, lty=3)
  # lines(x,Z, col='red')
  
  
  
  
  f    <- function(x,O,Sm,v) O[1]*x+O[2]*(1-exp(-O[3]*x))*Sm-v
  dfdx <- function(x,O,Sm) O[1]+O[2]*O[3]*exp(-O[3]*x)*Sm 
  ye   <- function(O, v){
    
    n  = length(v)
    w  = numeric(n)
    Sm = 1
    
    i=1
    w[i] = v[i]/O[1]
    for(i in 2:n){
      x0 = 0
      xdiff = 1
      tol = 1e-6
      iter = 0
      while(xdiff > tol & iter < 10){
        x1 = x0-f(x0,O,Sm,v[i])/dfdx(x0,O,Sm)
        xdiff = x1-x0
        x0 = x1
        iter = iter+1
        # cat("\niter",iter,"x0:",x0)
      }
      # w[i] = uniroot(f, c(0,100))$root
      w[i] = x0
      Sm = exp(-O[3]*w[i])*Sm+1
    }
    
    out = cumsum(w)
    return(out)
  }
  
  
  set.seed(1234)
  n = 100
  v  = rexp(n) # sample events in S-time
  pars = c(0.018, 0.74, 0.021)
  
  data.frame(rib = rib(pars,cumsum(v), dt=.5), ye(pars, v), exp_rib =  exp_rib(pars,v))
  
  
  
*/
  