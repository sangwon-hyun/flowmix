#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat build_mn_long_rcpp(const arma::cube& mu,
                             const Rcpp::IntegerVector& ntlist,
                             int iclust_r) {

  int iclust = iclust_r - 1; // Convert from R 1-indexing to C++ 0-indexing
  int TT = mu.n_rows;        // Assuming mu is (TT x dimdat x numclust)
  int dimdat = mu.n_cols;
  int total_n = sum(ntlist);

  // Allocate the final matrix size exactly once
  arma::mat mn_long(total_n, dimdat);

  int row_offset = 0;
  for(int t = 0; t < TT; ++t) {
    int n_t = ntlist[t];
    // Grab the mean vector for this cluster and time point
    arma::rowvec current_mu = mu.slice(iclust).row(t);

    // Fill the block for this time point using Armadillo's efficient fill
    // This is much faster than R's rep() + rbind()
    mn_long.rows(row_offset, row_offset + n_t - 1).each_row() = current_mu;

    row_offset += n_t;
  }

  return mn_long;
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List prepare_vectorized_data_rcpp(const Rcpp::List& ylist,
                                        const arma::cube& mu,
                                        int iclust_r) {

  int TT = ylist.size();
  int dimdat = mu.n_cols;
  int iclust = iclust_r - 1; // R to C++ indexing

  // 1. Calculate total rows (N)
  int total_n = 0;
  std::vector<int> n_t_vec(TT);
  for(int t = 0; t < TT; ++t) {
    arma::mat temp = Rcpp::as<arma::mat>(ylist[t]);
    n_t_vec[t] = temp.n_rows;
    total_n += n_t_vec[t];
  }

  // 2. Pre-allocate both giant matrices
  arma::mat y_long(total_n, dimdat);
  arma::mat mn_long(total_n, dimdat);

  // 3. Fill matrices in a single pass
  int row_offset = 0;
  for(int t = 0; t < TT; ++t) {
    int n_t = n_t_vec[t];
    if(n_t > 0) {
      // Fill Y: Convert R matrix to Armadillo and place it in the big block
      y_long.rows(row_offset, row_offset + n_t - 1) = Rcpp::as<arma::mat>(ylist[t]);

      // Fill MN: Grab the row and repeat it
      arma::rowvec current_mu = mu.slice(iclust).row(t);
      mn_long.rows(row_offset, row_offset + n_t - 1).each_row() = current_mu;

      row_offset += n_t;
    }
  }

  // Return as a list so R can pass them to dmvnorm_fast
  return Rcpp::List::create(
    Rcpp::Named("y_long") = y_long,
    Rcpp::Named("mn_long") = mn_long
  );
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::NumericVector dmvnorm_long_fused_cpp(const arma::mat& y_long,
                                           const arma::cube& mu,
                                           const Rcpp::IntegerVector& ntlist,
                                           int iclust_r,
                                           const arma::mat& sigma) {

  int iclust = iclust_r - 1; // R to C++ indexing
  int total_n = y_long.n_rows;
  int TT = ntlist.size();
  int dimdat = y_long.n_cols;

  // 1. Pre-calculate the Cholesky/Inverse once for this cluster
  // This is the "hoisting" we discussed
  arma::mat rooti = arma::inv(arma::trimatu(arma::chol(sigma)));
  double rootisum = arma::sum(arma::log(rooti.diag()));
  double const log2pi = std::log(2.0 * M_PI);
  double const other_terms = rootisum - (double)dimdat / 2.0 * log2pi;

  Rcpp::NumericVector out(total_n);
  int row_idx = 0;

  // 2. The Fused Loop
  for(int t = 0; t < TT; ++t) {
    int n_t = ntlist[t];
    if(n_t == 0) continue;

    // Grab mean vector for this time point
    arma::rowvec mn = mu.slice(iclust).row(t);

    for(int i = 0; i < n_t; ++i) {
      // Direct math on the row to stay in L2 Cache
      // We subtract the mean and rotate in one temporary vector
      arma::rowvec z = (y_long.row(row_idx) - mn) * rooti;

      double quad = 0;
      for(int j = 0; j < 3; ++j) { // Hardcoded 3 for speed
        quad += z[j] * z[j];
      }

      out[row_idx] = std::exp(other_terms - 0.5 * quad);
      row_idx++;
    }
  }

  return out;
}

// // [[Rcpp::export]]
// Rcpp::List calculate_all_densities_tiled(SEXP y_field_ptr,
//                                          const arma::cube& mu,
//                                          const Rcpp::List& sigma_list,
//                                          int numclust,
//                                          int TT) {

//   // Get pointer to the field of matrices
//   Rcpp::XPtr<arma::field<arma::mat>> ptr(y_field_ptr);
//   arma::field<arma::mat>& y_field = *ptr;

//   Rcpp::List final_denslist(numclust);

//   // Inside your function...
//   int dimdat = y_long.n_cols; // or 3 if you are hardcoding
//   double const log2pi = std::log(2.0 * M_PI); // Define it here

//   for (int iclust = 0; iclust < numclust; ++iclust) {
//     Rcpp::List dens_at_time(TT);

//     // 1. Hoist cluster-level linear algebra
//     arma::mat sigma = Rcpp::as<arma::mat>(sigma_list[iclust]);
//     arma::mat const rooti = arma::inv(arma::trimatu(arma::chol(sigma)));
//     double const rootisum = arma::sum(arma::log(rooti.diag()));
//     double const constants = -(double)3.0 / 2.0 * log2pi; // dimdat = 3
//     double const other_terms = rootisum + constants;

//     for (int t = 0; t < TT; ++t) {
//       const arma::mat& Y = y_field(t);
//       int n = Y.n_rows;
//       Rcpp::NumericVector out(n);
//       arma::rowvec mn = mu.slice(iclust).row(t);

//       // // 2. The Math Core using your optimized helper
//       // for (int i = 0; i < n; i++) {
//       //   arma::rowvec z = (Y.row(i) - mn);
//       //   inplace_tri_mat_mult(z, rooti);
//       //   // We stay in log-space as long as possible for stability
//       //   out[i] = std::exp(other_terms - 0.5 * arma::dot(z, z));
//       // }

//       // Try replacing the inner loop math with this:
//       for (int i = 0; i < n; i++) {
// 	// Armadillo's internal 'trimat' logic is quite fast
// 	// and might trigger the vectorizer better than a manual loop
// 	arma::rowvec z = (Y.row(i) - mn) * rooti;
// 	out[i] = std::exp(other_terms - 0.5 * arma::dot(z, z));
//       }

//       dens_at_time[t] = out;
//     }
//     final_denslist[iclust] = dens_at_time;
//   }
//   return final_denslist;
// }
