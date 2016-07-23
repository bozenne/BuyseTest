#include <iostream>
#include <RcppArmadillo.h>
#include <Rmath.h>
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp ;
using namespace std ;
using namespace arma ;

// [[Rcpp::export]]
mat createMat(){
mat A = randu<mat>(5,10);

vector<int> rm(3);
rm[0] = 1;
rm[1] = 2;
rm[3] = 3;

//A.shed_row(2);
//A.shed_cols(2,4);
A.shed_cols(rm);

return(A);
}