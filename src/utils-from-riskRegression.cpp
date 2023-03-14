// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

//' Row-wise cumulative sum 
//'
//' @description Fast computation of apply(x,1,cumsum)
//' @param X A matrix.
//' @return A matrix of same size as x.
// [[Rcpp::export(".rowCumSum_cpp")]]
arma::mat rowCumSum_cpp(const arma::mat X){
  arma::mat result = cumsum(X,1);
  return result;
}

//' Column-wise cumulative sum 
//'
//' @description Fast computation of apply(x,2,cumsum)
//' @param X A matrix.
//' @return A matrix of same size as x.
// [[Rcpp::export(".colCumSum_cpp")]]
arma::mat colCumSum_cpp(const arma::mat X){
  arma::mat result = cumsum(X,2);
  return result;
}

//' Apply cumprod in each row 
//'
//' @description Fast computation of t(apply(x,1,cumprod))
//' @param X A matrix.
//' @return A matrix of same size as x.
// [[Rcpp::export(".rowCumProd_cpp")]]
arma::mat rowCumProd_cpp(const arma::mat X){
  arma::mat result = cumprod(X,1);
  return result;
}

//' Substract a vector of values in each column 
//'
//' @description Fast computation of sweep(X, FUN = "-", STATS = center, MARGIN = 1)
//' @param X A matrix.
//' @param center A vector with length the number of rows of X .
//' @return A matrix of same size as x.
// [[Rcpp::export(".colCenter_cpp")]]
arma::mat colCenter_cpp(const arma::mat X, const arma::colvec& center){
  arma::mat result = X;
  result.each_col() -= center;
  return(result);
}

//' Substract a vector of values in each row
//'
//' @description Fast computation of sweep(X, FUN = "-", STATS = center, MARGIN = 2)
//' @param X A matrix.
//' @param center A vector with length the number of columns of X.
//' @return A matrix of same size as x.
// [[Rcpp::export(".rowCenter_cpp")]]
arma::mat rowCenter_cpp(const arma::mat X, const arma::rowvec& center){
  arma::mat result = X;
  result.each_row() -= center;
  return(result);
}

//' Divide by a vector of values in each column 
//'
//' @description Fast computation of sweep(X, FUN = "/", STATS = scale, MARGIN = 1)
//' @param X A matrix.
//' @param scale A vector with length the number of rows of X .
//' @return A matrix of same size as x.
// [[Rcpp::export(".colScale_cpp")]]
arma::mat colScale_cpp(const arma::mat X, const arma::colvec& scale){
  arma::mat result = X;
  result.each_col() /= scale;
  return(result);
}

//' Dividy by a vector of values in each row
//'
//' @description Fast computation of sweep(X, FUN = "/", STATS = center, MARGIN = 2)
//' @param X A matrix.
//' @param scale A vector with length the number of columns of X.
//' @return A matrix of same size as x.
// [[Rcpp::export(".rowScale_cpp")]]
arma::mat rowScale_cpp(const arma::mat X, const arma::rowvec& scale){
  arma::mat result = X;
  result.each_row() /= scale;
  return(result);
}

//' Multiply by a vector of values in each column 
//'
//' @description Fast computation of sweep(X, FUN = "*", STATS = scale, MARGIN = 1)
//' @param X A matrix.
//' @param scale A vector with length the number of rows of X .
//' @return A matrix of same size as x.
// [[Rcpp::export(".colMultiply_cpp")]]
arma::mat colMultiply_cpp(const arma::mat X, const arma::colvec& scale){
  arma::mat result = X;
  result.each_col() %= scale;
  return(result);
}

//' Multiply by a vector of values in each row
//'
//' @description Fast computation of sweep(X, FUN = "*", STATS = center, MARGIN = 2)
//' @param X A matrix.
//' @param scale A vector with length the number of columns of X.
//' @return A matrix of same size as x.
// [[Rcpp::export(".rowMultiply_cpp")]]
arma::mat rowMultiply_cpp(const arma::mat X, const arma::rowvec& scale){
  arma::mat result = X;
  result.each_row() %= scale;
  return(result);
}
