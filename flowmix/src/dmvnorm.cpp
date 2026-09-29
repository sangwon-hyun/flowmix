// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>

static double const log2pi = std::log(2.0 * M_PI);

// BRUTE FORCE: Tell the compiler to use O3 even if R passed O2
#pragma GCC optimize ("O3")

/* C++ version of the dtrmv BLAS function */
void inplace_tri_mat_mult(arma::rowvec &x, arma::mat const &trimat){
  arma::uword const n = trimat.n_cols;

  for(unsigned j = n; j-- > 0;){
    double tmp(0.);
    for(unsigned i = 0; i <= j; ++i)
      tmp += trimat.at(i, j) * x[i];
    x[j] = tmp;
  }
}

// [[Rcpp::export]]
arma::vec dmvnorm_arma_fast(arma::mat const &x,
			    arma::rowvec const &mean,
			    arma::mat const &sigma,
			    bool const logd = false) {
    using arma::uword;
    uword const n = x.n_rows,
             xdim = x.n_cols;
    arma::vec out(n);

    arma::mat const rooti = arma::inv(trimatu(arma::chol(sigma)));
    double const rootisum = arma::sum(log(rooti.diag())),
                constants = -(double)xdim/2.0 * log2pi,
              other_terms = rootisum + constants;

    arma::rowvec z;
    for (uword i = 0; i < n; i++) {
        z = (x.row(i) - mean);
        inplace_tri_mat_mult(z, rooti);
        out(i) = other_terms - 0.5 * arma::dot(z, z);
    }

    if (logd)
      return out;
    return exp(out);
}





// All credit goes to https://gallery.rcpp.org/articles/dmvnorm_arma/


// This is just a version that takes in rooti and other_terms, having done the
// cholesky decomposition of the covariance matrix sigma.
//[[Rcpp::export]]
arma::vec dmvnorm_arma_fast_kernel(arma::mat const &x,
			arma::rowvec const &mean,
			arma::mat const &rooti,
			double const other_terms) {

// This assumes xdim = 3

  using arma::uword;
  uword const n = x.n_rows;
  arma::vec out(n);

  // Use a standard C++ array/vector for the small, temporary result vector z
  arma::vec z(3);
  double density_contribution;
  double z_dot_z;

  for (uword i = 0; i < n; i++) {
      z_dot_z = 0.0;

      // 1. Subtraction and Multiplication (z = (x.row(i) - mean) * rooti)
      // Perform the full operation (z' = (x - mu) * R) using simple C++ loops
      // where R is rooti (D x D) and (x - mu) is 1 x D.

      for (int j = 0; j < 3; j++) { // Loop over dimensions (columns of rooti)
          // Calculate the j-th component of the resulting vector z'
          // z'[j] = sum_k (x[i,k] - mean[k]) * rooti[k, j]

          density_contribution = 0.0;
          for (int k = 0; k < 3; k++) { // Inner loop for dot product
              density_contribution += (x.at(i, k) - mean.at(k)) * rooti.at(k, j);
          }
          z(j) = density_contribution;
          z_dot_z += z(j) * z(j); // Calculate the squared magnitude (z.dot(z))
      }

      // 2. Final Density Calculation
      out(i) = other_terms - 0.5 * z_dot_z;
  }

  return exp(out); // or return out; if logd is true
}





// This is an older, more readable, but slightly slower, version of the above.

// // [[Rcpp::export]]
// arma::vec dmvnorm_arma_fast_kernel(arma::mat const &x,
// 			           arma::rowvec const &mean,
// 				   arma::mat const &rooti,
// 				   double const other_terms){
//     using arma::uword;
//     uword const n = x.n_rows,
//              xdim = x.n_cols;
//     arma::vec out(n);


//     arma::rowvec z;
//     for (uword i = 0; i < n; i++) {
//         z = (x.row(i) - mean);
//         inplace_tri_mat_mult(z, rooti);
//         out(i) = other_terms - 0.5 * arma::dot(z, z);
//     }

//     return exp(out);
// }
