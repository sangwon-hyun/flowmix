#include <numeric>
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat estepC(const arma::mat& ylong,
		 const arma::mat& mnlong,
		 const arma::vec& sqrt_resp_long,
		 const double& resp_sum
		 ){

  arma::mat resid = ylong - mnlong;
  resid.each_col() %= sqrt_resp_long;
  return ( resid.t() * resid ) / resp_sum;
  // return resid;

}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat ss(const arma::mat& X,
	     const arma::uvec& ind) {

              // int n = X.n_rows, k = X.n_cols;
              // arma::mat X(X.begin(), n, k, false);
              // arma::uvec ind = as<arma::uvec>(ind_);
              arma::mat submat = X.rows(ind - 1);
              // return wrap(submat);
	      return submat;
}
