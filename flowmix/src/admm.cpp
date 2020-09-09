#include <numeric>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector soft_threshC(NumericVector a,  double b) {
  LogicalVector mysign = (a > 0);
  return (as<NumericVector>(mysign) - 0.5) * 2 * pmax(abs(a) - b, 0);
}


// [[Rcpp::export]]
NumericVector wvec_updateC(NumericVector b1, NumericVector uw,
		    double lambda, double rho) {
  return soft_threshC(b1 + uw/rho, lambda/rho);
}


// [[Rcpp::export]]
NumericVector rowSumsC(NumericMatrix x) {
  int nrow = x.nrow(), ncol = x.ncol();
  NumericVector out(nrow);

  for (int i = 0; i < nrow; i++) {
    double total = 0;
    for (int j = 0; j < ncol; j++) {
      total += x(i, j);
    }
    out[i] = total;
  }
  return out;
}


// [[Rcpp::export]]
NumericVector rowSumsC_arma(arma::mat& x) {
  int nrow = x.n_rows, ncol = x.n_cols;
  NumericVector out(nrow);

  for (int i = 0; i < nrow; i++) {
    double total = 0;
    for (int j = 0; j < ncol; j++) {
      total += x(i, j);
    }
    out[i] = total;
  }
  return out;
}


// [[Rcpp::export]]
arma::mat projCmatC(arma::mat mat,
		    double C){

  arma::mat mat2 = mat % mat;
  NumericVector vlens = sqrt(rowSumsC_arma(mat2));
  vlens = pmax(vlens, C) / C;
  int nrow = mat.n_rows, ncol = mat.n_cols;
  arma::mat out(nrow, ncol);
  for(int ii = 0; ii < nrow; ii++) {
    out.row(ii) = mat.row(ii) / vlens[ii];
  }
  return out;
}


// [[Rcpp::export]]
arma::mat Z_updateC(arma::mat Xbeta1,
		    arma::mat Uz,
		    double C,
		    double rho,
		    int dimdat,
		    int TT) {
  return projCmatC(Xbeta1 + Uz / rho, C);
}
