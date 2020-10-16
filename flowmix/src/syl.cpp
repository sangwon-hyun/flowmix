#include <numeric>
#include <RcppArmadillo.h>
#include <RcppEigen.h>
using namespace Rcpp;
using namespace arma;
using namespace Eigen;

using Eigen::Map;                       // 'maps' rather than copies
using Eigen::MatrixXd;                  // variable size matrix, double precision

// [[Rcpp::export]]
Eigen::MatrixXd matrix_function_solve_triangular_sylvester_barebones(const Eigen::MatrixXd & TA,
								     const Eigen::MatrixXd & TB,
								     const Eigen::MatrixXd & C){
  // Eigen::eigen_assert(TA.rows() == TA.cols());
  // Eigen::eigen_assert(TA.Eigen::isUpperTriangular());
  // Eigen::eigen_assert(TB.rows() == TB.cols());
  // Eigen::eigen_assert(TB.Eigen::isUpperTriangular());
  // Eigen::eigen_assert(C.rows() == TA.rows());
  // Eigen::eigen_assert(C.cols() == TB.rows());

  // typedef typename MatrixType::Index Index;
  // typedef typename MatrixType::Scalar Scalar;

  int m = TA.rows();
  int n = TB.rows();
  Eigen::MatrixXd X(m, n);

  for (int i = m - 1; i >= 0; --i) {
    for (int j = 0; j < n; ++j) {

      // Compute T_A X = \sum_{k=i+1}^m T_A_{ik} X_{kj}
      double TAX;
      if (i == m - 1) {
      	TAX = 0;
      } else {
	MatrixXd TAXmatrix = TA.row(i).tail(m-1-i) * X.col(j).tail(m-1-i);
      	TAX = TAXmatrix(0,0);
      }

      // Compute X T_B = \sum_{k=1}^{j-1} X_{ik} T_B_{kj}
      double XTB;
      if (j == 0) {
      	XTB = 0;
      } else {
      	MatrixXd XTBmatrix = X.row(i).head(j) * TB.col(j).head(j);
      	XTB = XTBmatrix(0,0);
      }

      X(i,j) = (C(i,j) - TAX - XTB) / (TA(i,i) + TB(j,j));
    }
  }
  return X;
}



// // [[Rcpp::export]]
// Eigen::MatrixXd sylC_upper_tri_barebones(const Eigen::MatrixXd & TA,
// 					 const Eigen::MatrixXd & TB,
// 					 const Eigen::MatrixXd & F){
//   Eigen::MatrixXd Y = Eigen::internal::matrix_function_solve_triangular_sylvester(TA, TB, F);
// }

// [[Rcpp::export]]
Eigen::MatrixXd sylC_upper_tri(const Eigen::MatrixXd & TA,
			       const Eigen::MatrixXd & UA,
			       const Eigen::MatrixXd & tUA,
			       const Eigen::MatrixXd & TB,
			       const Eigen::MatrixXd & UB,
			       const Eigen::MatrixXd & tUB,
			       const Eigen::MatrixXd & syl_C){

  // Form the constant
  Eigen::MatrixXd F = ((tUA * syl_C) * UB).real();
  return F;
  // Eigen::MatrixXd F = (-1) * (UA.adjoint() * syl_C) * UB;

  // // // Solve the sylvester equation given upper triangular matrices R and S
  // Eigen::MatrixXd Y = Eigen::internal::matrix_function_solve_triangular_sylvester(TA, TB, F);

  // // // convert back to original scale
  // Eigen::MatrixXd X = (UA * Y) * tUB();
  // // // // Eigen::MatrixXd X = (UA * Y) * UB.transpose();

  // return X.transpose();
}


// [[Rcpp::export]]
arma::mat prepare_sylC_const3(const arma::mat& U,
			      const arma::mat& Xaug,
			      const double& rho,
			      const arma::mat& Z,
			      const arma::mat& X,
			      const arma::mat& W,
			      const arma::mat& term3,
			      const arma::mat& sigma,
			      const arma::mat& Xinv){

			      // const int N,
			      // const arma::mat& sigmainv,
			      // const arma::mat& yXcentered,

  // Setup (SH: correct)
  arma::mat term1 = U.t() * Xaug;//t(U) %*% Xaug;
  arma::mat term2 = rho * (Z.t() * X + W.t());//rho *  (t(Z) %*% X + t(W));
  // arma::mat term3 = sigmainv * yXcentered;//1/N * sigmainv %*% yXcentered;
  // term3 = term3 / N;
  return((sigma * (term1 - term2 - term3 ) * Xinv));
}


// [[Rcpp::export]]
arma::mat prepare_sylC_const33(const arma::mat& U,
			      const arma::mat& Xaug,
			      const double& rho,
			      const arma::mat& Z,
			      const arma::mat& X,
			      const arma::mat& W,
			      const int& N,
			      const arma::mat& sigmainv,
			      const arma::mat& yXcentered,
			      const arma::mat& sigma,
			      const arma::mat& Xinv){

			      // const int N,
			      // const arma::mat& sigmainv,
			      // const arma::mat& yXcentered,

  // Setup (SH: correct)
  arma::mat term1 = U.t() * Xaug;//t(U) %*% Xaug;
  arma::mat term2 = rho * (Z.t() * X + W.t());//rho *  (t(Z) %*% X + t(W));
  arma::mat term3 = sigmainv * yXcentered;//1/N * sigmainv %*% yXcentered;
  // term3 = term3 / N;
  return( sigma * (term1 - term2 - term3 ) * Xinv);
}


// Taken from:
// https://discourse.mc-stan.org/t/solve-a-lyapunov-sylvester-equation-include-custom-c-function-using-eigen-library-possible/12688
// Also helpful https://stackoverflow.com/questions/56929966/implementing-the-bartels-stewart-algorithm-in-eigen3
