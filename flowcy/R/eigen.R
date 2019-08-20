## Synopsis: Do eigendecomposition once, use thrice, use thrice


##' From eigendecomposition of the sigmas, calculate dmvnorm.
##' @param sigma_eig Result of eigendecomposition of sigma.
dmvnorm_fast <- function(y, mu, sigma_eig){
  const = (2 * pi)^(- dimdat/2)

  ## Setup
  TT = nrow(y)
  mydet <- sigma_eig$det
  myinv <- sigma_eig$sigma_inv
  dimdat <- ncol(y)

  ## Main calculations.
  resids = sweep(y, 2, mu) ## This takes a lot of time, but not clear how to fix
  myinv_half <- sigma_eig$inverse_sigma_half
  transformed_resids = resids %*% myinv_half

  ## Twice as fast as regular RowSums(), which is supposed to be fast itself
  ## resid.prods = Rfast::rowsums(transformed_resids * transformed_resids)
  resid.prods = rowSums(transformed_resids * transformed_resids)

  return(const * mydet^(-1/2) * exp(-1/2 * resid.prods))
}



##' Get's eigendecomposition of a matrix. Basically, \code{sigma == evecs %*% Lambdamat
##' %*% solve(evecs)}.
##' @param sigma dimdat by dimdat matrix
##' @return List containing eigenvalues (vector), eigenvectors (matrix of
##'   eigenvectors as columns)
eigendecomp_sigma_barebones <- function(sigma){

  ## Do eigendecomposition
  eig = eigen(sigma)

  ## Gather and return results
  sigma_eig <- list(values=eig$values, vectors=eig$vectors)
  return(sigma_eig)
}



##' Get's eigendecomposition of a matrix. Basically, \code{sigma == evecs %*% Lambdamat
##' %*% solve(evecs)}.
##' @param sigma A single (dimdat x dimdat) matrix.
##' @return List containing eigenvalues (vector), eigenvectors (matrix of
##'   eigenvectors as columns)
eigendecomp_sigma <- function(sigma){

  ## Do a barebones eigendecomposition (this is the part that can be sped up
  ## using armadillo R)
  sigma_eig0 = eigendecomp_sigma_barebones(sigma)

  ## Gather and return results
  sigma_eig <- list(values = sigma_eig0$values,
                    vectors = sigma_eig0$vectors,
                    det = det_from_eig(sigma_eig0),
                    sigma_half = sigma_half_from_eig(sigma_eig0),
                    inverse_sigma_half = inverse_sigma_half_from_eig(sigma_eig0),
                    sigma_inv = sigma_inv_from_eig(sigma_eig0),
                    sigma = sigma
                    )
  return(sigma_eig)
}


##' From eigendecomposition of the sigmas, calculate determinant.
##' @param sigma_eig Eigendecomposition (from \code{get_sigma_eig()}) of sigma.
##' @return Determinant of sigma.
det_from_eig <- function(sigma_eig){
  return(prod(sigma_eig$values))
}


##' From the eigendecomposition, make inverse of sigma.
##' @param sigma_eig Eigendecomposition (from \code{get_sigma_eig()}) of sigma.
##' @return Inverse of sigma.
sigma_inv_from_eig <- function(sigma_eig){
  (sigma_eig$vectors %*% diag(1/sigma_eig$values) %*% t(sigma_eig$vectors))
}


##' From the eigendecomposition, make sigma halves.
##' @param sigma_eig Eigendecomposition (from \code{get_sigma_eig()}) of sigma.
##' @return Square root of sigma i.e. if the returned object is
##'   \code{sigma_half}, we should have \code{sigma_half %*% sigma_half ==
##'   sigma}.
sigma_half_from_eig <- function(sigma_eig){
  (sigma_eig$vectors %*% diag(sqrt(sigma_eig$values)) %*% t(sigma_eig$vectors))
}


##' From the eigendecomposition, make inverse sigma halves.
##' @param sigma_eig Eigendecomposition (from \code{get_sigma_eig()}) of sigma.
##' @return Square root of sigma i.e. if the returned object is
##'   \code{sigma_half}, we should have \code{sigma_half %*% sigma_half ==
##'   sigma}.
inverse_sigma_half_from_eig <- function(sigma_eig){
  (sigma_eig$vectors %*% diag(1/sqrt(sigma_eig$values)) %*% t(sigma_eig$vectors))
} ## TODO Test this.

##' From a (TT x numclust x dimdat x dimdat) array whose [tt,ii,,]'th entry is
##' the (dimdat x dimdat) covariance matrix, do eigendecompositions.
##' @param sigma_array (TT x numclust x dimdat x dimdat) array
##' @return TT length list of (numclust length list of eigendecompositions).
eigendecomp_sigma_array <- function(sigma_array){

  TT = dim(sigma_array)[1]
  numclust = dim(sigma_array)[2]

  ## Only need to calculate once because sigmas are the same across tt=1:TT
  eig_by_dim = lapply(1:numclust, function(idim){
    eigendecomp_sigma(sigma_array[1, idim, , ])
  })
  ## eig_list = lapply(1:TT, function(tt){ eig_by_dim }) ## In fact, it would be
  ##                                                     ## better if we just
  ##                                                     ## passed around 1 thing

  ## ## The old, incredibly repetitive way of doing it:
  ## for(iclust in 1:numclust){
  ##   eig_list[[tt]] = lapply(1:numclust, function(idim){
  ##     eigendecomp_sigma(sigma_array[tt, idim, , ])
  ##   })
  ## }

  ## Todo: eventually, figure out a way to save sigma in a non-repetitive way.

  return(eig_by_dim)
}
