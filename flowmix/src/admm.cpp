#include <numeric>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
NumericVector soft_threshC(const NumericVector& a,  const double& b) {
  LogicalVector mysign = (a > 0);
  return (as<NumericVector>(mysign) - 0.5) * 2 * pmax(abs(a) - b, 0);
}


// [[Rcpp::export]]
NumericVector wvec_updateC(const NumericVector& b1,
			   const NumericVector& uw,
			   const double& lambda,
			   const double& rho) {
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
NumericVector rowSumsC_arma(const arma::mat& x) {
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
NumericVector rowSumsC2_arma(const arma::mat& x) {
  arma::mat y = x % x;
  int nrow = y.n_rows, ncol = y.n_cols;
  NumericVector out(nrow);

  for (int i = 0; i < nrow; i++) {
    double total = 0;
    for (int j = 0; j < ncol; j++) {
      total += y(i, j);
    }
    out[i] = total;
  }
  return out;
}


// [[Rcpp::export]]
arma::mat projCmatC(const arma::mat& mat,
		    const double& C){

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




// [[Rcpp::export]]
arma::mat mv_mult(const arma::mat& lhs,
                  const arma::vec& rhs)
{
  return lhs * rhs;
}


//[[Rcpp::export]]
arma::vec b_updateC(const NumericVector& wvec,
		    const NumericVector& uw,
		    const double& rho,
		    const NumericVector& cvec3_el,
		    const NumericVector& yvec,
		    const arma::mat& DtDinv,
		    const arma::mat& D,
		    const double& N){

  // Calculate the three subvectors
  arma::vec cvec1 = sqrt(1/(2*N)) * yvec;
  arma::vec cvec2 = sqrt(rho/2) * (wvec - uw/rho);
  arma::vec cvec3 = sqrt(rho/2) * cvec3_el;
  arma::vec cvec;
  cvec = arma::join_cols(cvec1, cvec2, cvec3);

  // Do the solve.
  arma::mat betahat = DtDinv * D.t() * cvec;
  return betahat;
}


//[[Rcpp::export]]
arma::mat subtractC2(const arma::vec& wt,
		     const arma::mat& mat,
		     const arma::vec& vec,
		     const arma::mat& mat2
		     ){
  int nr = mat.n_rows;
  int nc = mat.n_cols;
  arma::mat vect = vec.t();
  arma::mat matnew(nr, nc);
  arma::rowvec rowvec(nc);
  for(int ii = 0; ii < nr; ii++){
    matnew.row(ii) = wt(ii) * (mat.row(ii) - vect) * mat2;
  }
  return matnew;
}



//[[Rcpp::export]]
arma::mat subtractC3(const arma::vec& wt,
		     const arma::mat& mat,
		     const arma::vec& vec
		     ){
  int nr = mat.n_rows;
  int nc = mat.n_cols;
  arma::mat vect = vec.t();
  arma::mat matnew(nr, nc);
  arma::rowvec rowvec(nc);
  for(int ii = 0; ii < nr; ii++){
    matnew.row(ii) = wt(ii) * (mat.row(ii) - vect);
  }
  return matnew;
}


//[[Rcpp::export]]
arma::mat dothisC(const arma::vec& longwt,
		  const arma::mat& ylong,
		  const arma::mat& mumat,
		  const arma::mat& sigma_half){

  arma::mat matnew = (ylong - mumat) * sigma_half;
  int nr = ylong.n_rows;
  for(int ii = 0; ii < nr; ii++){
    matnew.row(ii) = longwt(ii) * matnew.row(ii);
  }
  return rowSumsC2_arma(matnew);
}
// ## wtresidmat = longwt * (ylong - mumat)
//  ## transformed_resids = wtresidmat %*% sigma_half
//   ## resid.prods = rowSumsC2_arma(transformed_resids)
    // resid.prods = dothisC(longwt, ylong, mumat, sigma_half)
