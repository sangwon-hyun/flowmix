##' (Helper) Soft thresholding of |a| at radius of |b|.
soft_thresh <- function(a, b){
  return(sign(a) * pmax(0, abs(a)-b))

  ## (Kind of) test for this function:
  x = seq(from=-10, to=10, length=100)
  plot(y=soft_thresh(x, 3), x=x)
  abline(h=0)
  abline(v=c(-3,3), col='grey')
  abline(v=0)
}


##' (Helper) Projection of each ROW of the matrix \code{mat} into an l2 ball
##' centered at zero, of radius C.
projCmat <- function(mat, C){
  if(!is.null(C)){
    vlens = sqrt(rowSums(mat * mat))
    inds = which(vlens > C)
    if(length(inds) > 0){
      mat[inds,] = mat[inds,] * C / vlens[inds]
    }
  }
  return(mat)
}

##' (Helper) Projection into an l2 ball centered at zero, of radius C.
projC <- function(v, C){
  vlen = sqrt(sum(v*v))
 if(vlen <  C){
    return(v)
  } else {
    return(C * v/vlen)
  }

  ## ## Test for this function:
  ## maxdev = 0.5
  ## for(isim in 1:1000){
  ##   set.seed(isim)
  ##   v = rnorm(10)
  ##   vlen = sqrt(sum(v*v))
  ##   vproj = projC(v, maxdev)
  ##   if(vlen >= maxdev){
  ##     assert_that(all.equal(sqrt(sum(vproj * vproj)),maxdev)==TRUE)
  ##   }
  ##   assert_that(all(abs((vproj/v - (vproj/v)[1])) < 1E-15 ))
  ## }
}



##' (Helper) Check convergence for ADMM.
##' @param err_rel = 1E-3
converge <- function(beta1, rho, w, Z, w_prev, Z_prev, Uw, Uz, tX, Xbeta1,
                     err_rel = 1E-4
                     ## , err_abs = 1E-3
                     ){

  ## Form primal and dual residuals, and other things.
  ## prim1 = A %*% beta1 ## old
  ## prim2 = B %*% wz    ## old
  prim1 = rbind(beta1, Xbeta1) ## new (although I don't like the rbind)
  prim2 = - rbind(w, Z) ## new
  primal_resid = prim1 + prim2
  ## dual_resid = rho * tAB %*% (wz - wz_prev) ## old (No need for tAB).
  dual_resid = -rho * ((w - w_prev) + (tX %*% (Z - Z_prev))) ## new
  tAU = Uw + tX %*% Uz

  ## Form primal and dual tolerances.
  primal_err = ## sqrt(length(primal_resid)) * err_abs +
    err_rel * max(norm(prim1, "F"), norm(prim2, "F"))
  dual_err = ## sqrt(length(dual_resid)) * err_abs +
    err_rel * norm(tAU, "F")

  ## Check convergence.
  primal_converge = ( norm(primal_resid, "F") <= primal_err )
  dual_converge = ( norm(dual_resid, "F") <= dual_err )

  ## return(primal_converge & dual_converge)
  converge = primal_converge & dual_converge
  return(list(primal_resid = primal_resid,
              primal_err = primal_err,
              dual_resid = dual_resid,
              dual_err = dual_err,
              converge = converge))
}


##' (Helper) calculates the per-cluster objective value for the ADMM. INCOMPLETE
##' because it is super inefficient right now.
##' @param beta p x d matrix
objective_per_cluster <- function(beta, ylist, Xa, resp, lambda, dimdat, iclust, sigma, iter,
                                  first_iter, sigma_eig_by_clust=NULL, mc.cores = 1){

  ## Setup
  TT = length(ylist)
  p = ncol(Xa) - 1
  mn = Xa %*% beta

  ## Obtain square root of sigma.
  if(first_iter){
    sigma_half <- sigma_half_from_eig(eigendecomp_sigma_barebones(sigma))
  } else {
    sigma_eig <- sigma_eig_by_clust[[iclust]]
    sigma_half <- sigma_eig$inverse_sigma_half
  }

  ## Obtain the weighted residuals once.
  wt_resids_list = mclapply(1:TT, function(tt){
    y = ylist[[tt]]
    mumat = matrix(mn[tt,],
                   ncol = ncol(y),
                   nrow = nrow(y),
                   byrow = TRUE)
    wt_resid = sqrt(resp[[tt]][, iclust]) * (y - mumat)
  })

  ## Calculate the inner product resid * sigma^-1 * t(resid)
  transformed_resids = do.call(rbind, wt_resids_list) %*% sigma_half
  resid.prods = rowSums(transformed_resids * transformed_resids)
  grand.sum = sum(resid.prods)

  ## Calculate the objective value
  tol = 1E-8
  obj = (1/2) * grand.sum + lambda * sum(abs(beta) > tol)
  return(obj)
}


##' Calculate a specific least squares problem \min_b \|c-Db\|^2.
b_update  <- function(wvec, uw, Z, Uz, rho, yvec, D, DtDinv){

  cvec = c(sqrt(1/2) * yvec,
           sqrt(rho/2) * (wvec - uw/rho),
           sqrt(rho/2) * as.numeric(t(Z - Uz/rho)))
  sol = DtDinv %*% crossprod(D, cvec)
  return(sol)
}

##' Makes a permutation matrix with columns corresponding to |intercept_inds|
##' are all zero, and the rest are an identity matrix. left-multiplication of
##' this matrix times a (p+1)*dimdat vector like |b| extracts the non-intercept
##' entries.
##' @param p Number of covariates.
##' @param dimdat Data dimension.
##' @param intercept_inds Numeric vector containing indices of the intercept
##'   entries.
##'
##' (TODO: improve; this is super inefficient)
##' @return ((p dimdat) x ((p+1) dimdat)) permutation matrix.
make_I_aug <- function(p, dimdat, intercept_inds){
  I = diag(rep(1, p * dimdat))
  Imaster = matrix(0,
                   nrow = p * dimdat,
                   ncol = (p+1) * dimdat)
  Imaster[,-intercept_inds] = I
  return(Imaster)
}

  ## ## (Not used for now) Define a fast linear system solver
  ## Rcpp::cppFunction(depends='RcppArmadillo', code='
  ##         arma::mat fRcpp (arma::mat A, arma::mat b) {
  ##         arma::mat betahat ;
  ##         betahat = (A.t() * A ).i() * A.t() * b ;
  ##         return(betahat) ;
  ##         }
  ##         ')
