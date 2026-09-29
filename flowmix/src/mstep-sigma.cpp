#include <numeric>
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat mstep_sigma_C(const arma::mat& ylong,
		 const arma::mat& mnlong,
		 const arma::vec& sqrt_resp_long,
		 const double& resp_sum
		 ){

  arma::mat resid = (ylong - mnlong);

  // We should stop execution if the dimensions are wrong.
  if (sqrt_resp_long.n_rows != resid.n_rows) {
      Rcpp::stop("sqrt_resp_long must have the same number of rows as ylong.");
  }

  resid.each_col() %= sqrt_resp_long;
  return ( resid.t() * resid ) / resp_sum;
  // return resid;

}


