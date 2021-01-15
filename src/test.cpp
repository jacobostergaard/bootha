// #include <Rcpp.h>
// using namespace Rcpp;
#include <RcppArmadillo.h>
#include <stdio.h>
#include <math.h>

using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]
using namespace std;
using namespace arma;


// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
arma::vec testfun(arma::vec x){
  
  vec out = regspace(x(0),x(1),x(2));
    
  return out;
}


// ---------------- R code ---------------- //

/*** R

  # misc::clean_up()
  
  x = c(0,.2,10)
  testfun(x)
  
*/
