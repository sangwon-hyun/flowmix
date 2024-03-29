##' From eigendecomposition of the sigmas, calculate the same thing as
##' \code{mvtnorm::dmvnorm()}.
##'
##' @param y Multivariate data.
##' @param mu Mean vector.
##' @param sigma_eig Result of eigendecomposition of sigma (using
##'   \code{eigendecomp_sigma(sigma)} for covariance matrix \code{sigma}).
##'
##' @return Density vector.
dmvnorm_fast <- function(y, mu, sigma_eig){

  ## Basic checks
  assertthat::assert_that("matrix" %in% class(y))

  dimdat = ncol(y)
  const = (2 * pi)^(-dimdat/2)

  ## Setup
  TT = nrow(y)
  mydet <- sigma_eig$det
  myinv <- sigma_eig$sigma_inv
  dimdat <- ncol(y)

  ## Main calculations.
  mumat = matrix(mu,
                 ncol=ncol(y),
                 nrow=nrow(y),
                 byrow=TRUE)
  resids = y - mumat

  myinv_half <- sigma_eig$inverse_sigma_half
  transformed_resids = resids %*% myinv_half

  ## Twice as fast as regular RowSums(), which is supposed to be fast itself
  resid.prods = rowSums(transformed_resids * transformed_resids)
  dens = const * mydet^(-1/2) * exp(-1/2 * resid.prods)
  return(dens)
}



##' Gets eigendecomposition of a matrix. Basically, \code{sigma == evecs \%*\%
##' Lambdamat \%*\% solve(evecs)}.
##'
##' @param sigma dimdat by dimdat matrix
##'
##' @return List containing eigenvalues (vector), eigenvectors (matrix of
##'   eigenvectors as columns)
eigendecomp_sigma_barebones <- function(sigma){

  ## Do eigendecomposition
  eig = eigen(sigma)

  values = eig$values
  if(length(sigma)==1) values = as.matrix(values)

  ## Gather and return results
  sigma_eig <- list(values = values, vectors = eig$vectors)
  return(sigma_eig)
}



##' (Wrapper for \code{eigendecomp_sigma_barebones()}) Gets eigendecomposition
##' of a matrix. Basically, \code{sigma == evecs \%*\% Lambdamat \%*\%
##' solve(evecs)}.
##'
##' @param sigma A single (dimdat x dimdat) matrix.
##'
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
                    sigma = sigma)
  return(sigma_eig)
}


##' Gets cholesky decomposition.
##'
##' @param sigma A single (dimdat x dimdat) matrix.
##'
##' @return List containing eigenvalues (vector), eigenvectors (matrix of
##'   eigenvectors as columns)
choldecomp_sigma <- function(sigma){

  sigma_eig0 = eigendecomp_sigma_barebones(sigma)

  ## Gather and return results
  sigma_chol <- list(chol=chol(sigma),
                     values = sigma_eig0$values,
                     vectors = sigma_eig0$vectors,
                     det = det_from_eig(sigma_eig0),
                     sigma_half = sigma_half_from_eig(sigma_eig0),
                     inverse_sigma_half = inverse_sigma_half_from_eig(sigma_eig0),
                     sigma_inv = sigma_inv_from_eig(sigma_eig0),
                     sigma = sigma
                    )
}


##' From eigendecomposition of the sigmas, calculate determinant.
##' @param sigma_eig Eigendecomposition (from \code{get_sigma_eig()}) of sigma.
##' @return Determinant of sigma.
det_from_eig <- function(sigma_eig){
  return(prod(sigma_eig$values))
}


##' From the eigendecomposition, make inverse of sigma.
##'
##' @param sigma_eig Eigendecomposition (from \code{get_sigma_eig()}) of sigma.
##'
##' @return Inverse of sigma.
sigma_inv_from_eig <- function(sigma_eig){
  vec = 1 / sigma_eig$values
  if(length(vec)==1){
    mat = vec
  } else {
    mat = diag(vec)
  }
  (sigma_eig$vectors %*% mat %*% t(sigma_eig$vectors))
}


##' From the eigendecomposition, make sigma halves.
##'
##' @param sigma_eig Eigendecomposition (from \code{get_sigma_eig()}) of sigma.
##'
##' @return Square root of sigma i.e. if the returned object is
##'   \code{sigma_half}, we should have \code{sigma_half %*% sigma_half ==
##'   sigma}.
sigma_half_from_eig <- function(sigma_eig){
  vec = sqrt(sigma_eig$values)
  if(length(vec)==1){
    mat = vec
  } else {
    mat = diag(vec)
  }
  (sigma_eig$vectors %*% mat %*% t(sigma_eig$vectors))
}


##' From the eigendecomposition, make inverse sigma halves.
##'
##' @param sigma_eig Eigendecomposition (from \code{get_sigma_eig()}) of sigma.
##'
##' @return Square root of sigma i.e. if the returned object is
##'   \code{sigma_half}, we should have \code{all((sigma_half %*% sigma_half) ==
##'   sigma)}.
inverse_sigma_half_from_eig <- function(sigma_eig){
  vec = 1/sqrt(sigma_eig$values)
  if(length(vec)==1){
    mat = vec
  } else {
    mat = diag(vec)
  }
  (sigma_eig$vectors %*% mat %*% t(sigma_eig$vectors))
}

##' From a (numclust x dimdat x dimdat) array whose [ii,,]'th entry is
##' the (dimdat x dimdat) covariance matrix, do eigendecompositions.
##'
##' @param sigma_array (numclust x dimdat x dimdat) array.
##'
##' @return |numclust| lengthed list of eigendecompositions.
eigendecomp_sigma_array <- function(sigma_array){

  numclust = dim(sigma_array)[1]

  ## Only need to calculate once because sigmas are the same across tt=1:TT
  eig_by_clust = lapply(1:numclust, function(idim){
    eigendecomp_sigma(sigma_array[idim, , ])
  })

  return(eig_by_clust)
}


##' From a (TT x numclust x dimdat x dimdat) array whose [tt,ii,,]'th entry is
##' the (dimdat x dimdat) covariance matrix, do cholesky decompositions.
##' @param sigma_array (TT x numclust x dimdat x dimdat) array
##' @return TT length list of (numclust length list of eigendecompositions).
choldecomp_sigma_array <- function(sigma_array){

  TT = dim(sigma_array)[1]
  numclust = dim(sigma_array)[2]

  ## Only need to calculate once because sigmas are the same across tt=1:TT
  chol_by_clust = lapply(1:numclust, function(idim){
    choldecomp_sigma(sigma_array[1, idim, , ])
  })
  return(chol_by_clust)
}
