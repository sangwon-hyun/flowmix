#include <numeric>
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::mat mstep_sigma_indexed_C(const arma::mat& ylong,
				const arma::mat& mn_small,   // The T x P matrix (mn[,,iclust])
				const arma::uvec& ntlist,    // The number of points at each time T
				const arma::vec& sqrt_resp,  // The N-length vector of sqrt(weights)
				const double& resp_sum) {

  int dimdat = ylong.n_cols;
  arma::mat sigma = arma::zeros<arma::mat>(dimdat, dimdat);

  int current_start = 0;
  int TT = ntlist.n_elem;

  // Loop through each time point T
  for(int t = 0; t < TT; t++) {
    int n_t = ntlist(t);
    if(n_t == 0) continue;

    // 1. Get the block of data for this time point
    // 2. Subtract the single mean vector for this time point (mn_small.row(t))
    // 3. Weight by the corresponding block of sqrt_resp
    arma::mat resid_t = ylong.rows(current_start, current_start + n_t - 1);
    resid_t.each_row() -= mn_small.row(t);

    arma::vec weights_t = sqrt_resp.subvec(current_start, current_start + n_t - 1);
    resid_t.each_col() %= weights_t;

    // Accumulate the crossproduct
    sigma += resid_t.t() * resid_t;

    current_start += n_t;
  }

  return sigma / resp_sum;
}
