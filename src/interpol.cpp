// #include <Rcpp.h>
#include <RcppArmadillo.h>
#include <stdio.h>
#include <math.h>

using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]
using namespace std;
using namespace arma;

// Declare functions
arma::vec interpol_aux(arma::mat f, arma::vec x, bool flat=false) {
  
  int N = f.n_rows-1;     // number of input values in (x,f(x)) format to interpolate between
  int n = x.n_elem;       // number of values to interpolate at
  vec z = zeros<vec>(n);  // output vector of interpolated values
  
  uvec idx; // generic index vector
  double a; // generic slope value for non-flat extrapolation
  double b; // generic intercept value for non-flat extrapolation
  int id1;  // generic index value
  int id2;  // generic index value
  
  double x1; // generic x value
  double x2; // generic x value
  double y1; // generic y value
  double y2; // generic y value
  double w;  // weight for interpolation
  int j;     // generic index value
  
  
  idx = sort_index(f.col(0));  // make sure input values are sorted in x-value order
  f   = f.rows(idx);
  vec fx = f.col(0);
  vec fy = f.col(1);
  
  if(flat){ // extrapolate extremes flat
    idx    = find(x <= min(fx));
    z(idx) = ones<vec>(idx.n_elem)*fy(0);
    
    idx   = find(x >= max(fx));
    z(idx) = ones<vec>(idx.n_elem)*fy(N);
  } else{ // extrapolate extremes linearly
    
    idx    = find(x <= min(fx));
    a      = (fy(1)-fy(0))/(fx(1)-fx(0));
    b      = fy(1)-fx(1)*a;
    z(idx) = a*x(idx)+b;
    
    idx    = find(x >= max(fx));
    a      = (fy(N)-fy(N-1))/(fx(N)-fx(N-1));
    b      = fy(N)-fx(N)*a;
    z(idx) = a*x(idx)+b;
    
  }
  
  // interpolate linearly within extremes
  idx    = find(x > min(fx) && x < max(fx));
  
  for(int i=0;i<idx.n_elem;i++){ 
    j = idx(i);
    id1 = max(find( fx <= x(j) ));
    id2 = min(find( fx >= x(j) ));
    x1  = fx(id1);
    x2  = fx(id2);
    y1  = fy(id1);
    y2  = fy(id2);
    
    if(x2==x1){
      w = 1;
    } else{
      w = (x2-x(j))/(x2-x1);
    }
    z(j) = y1*w+y2*(1-w);
  }
  
  return(z);
}  // 


  // [[Rcpp::export]]
arma::vec interpol(arma::mat f, Rcpp::Nullable<NumericVector> x = R_NilValue, Rcpp::Nullable<NumericVector> y = R_NilValue, bool flat=false) {
  
  vec out;
  mat g = f;
  
  if (x.isNotNull()) {
    // interpolate from x to y
    out = interpol_aux(g,as<arma::vec>(x), flat);
    return(out);
    
  } else if(y.isNotNull()){
    // interpolate from y to x
    g.col(0) = f.col(1);
    g.col(1) = f.col(0);
    out = interpol_aux(g,as<arma::vec>(y), flat);
    
    return(out);
    
  } else{
    Rcout << "Either x or y values must be provided." << endl;
    return out;
  }
  
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//
  
  /*** R
    f = matrix(1:4,nr=2,nc=2)
    interpol(f,NULL,NULL)

    interpol(f,c(1,2,3),NULL)
    interpol(f,NULL,c(1,2,3))

    x = seq(1,10,.1)
    # y = sin(x)
    y = x^2+1
    f = cbind(x,y)

    plot(f, xlim=c(-10,20), ylim=c(-10,200), type='l', lwd=4)
    xin = c(-5,1,4,5,10,12)
    y1 = interpol(f,xin,NULL,TRUE)
    y2 = interpol(f,xin,NULL,FALSE)
    y3 = interpol(f,xin,NULL)

    points(xin,y1, pch=16, col='dodgerblue')
    points(xin,y2, pch=16, col='red')
    points(xin,y3, pch=16, col='orange')

    xin = seq(-10,20,.01)
    y4 = interpol(f,xin,NULL)
    lines(xin,y4, lty=1, col='red')



    yin = c(-5,10,50,100,150,175)
    interpol(f,y=yin,flat = TRUE)
    interpol(f,y=yin,flat = FALSE)

    (x1 = interpol(f,NULL,yin,TRUE))
    (x2 = interpol(f,NULL,yin,FALSE))
    (x3 = interpol(f,NULL,yin))



    points(x1,yin, pch=2, col='dodgerblue')
    points(x2,yin, pch=2, col='red')
    points(x3,yin, pch=2, col='orange')
  
    
*/
  