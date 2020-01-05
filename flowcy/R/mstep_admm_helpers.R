##' (Helper) Soft thresholding of |a| at radius of |b|.
soft_thresh <- function(a, b){
  return(sign(a) * pmax(0, abs(a)-b))

  ## ## (Kind of) test for this function:
  ## x = seq(from=-10, to=10, length=100)
  ## plot(y=soft_thresh(x, 3), x=x)
  ## abline(h=0)
  ## abline(v=c(-3,3), col='grey')
  ## abline(v=0)
}


##' (Helper) Projection of each ROW of the matrix \code{mat} into an l2 ball
##' centered at zero, of radius C.
projCmat <- function(mat, C){
  vlens = sqrt(rowSums(mat * mat))
  inds = which(vlens > C)
  if(length(inds)>0){
    mat[inds,] = mat[inds,] * C / vlens[inds]
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



##' (Helper) Check convergence for ADMM (incomplete because I couldn't figure
##' out what it should be when we are doing a matrix ADMM; will come back to
##' it).
##' @param err_rel = 1E-3
converge_old <- function(beta1, X, rho, w, Z, w_prev, Z_prev, Uw, Uz, err_rel = 1E-4,
                     A, B, tA, tAB
                     ## , err_abs = 1E-3
                     ){

  U = rbind(Uw, Uz)

  ## Form the second block of primal variables.
  wz = rbind(w,
             Z)
  wz_prev = rbind(w_prev,
                  Z_prev)

  ## Form primal and dual residuals.
  primal_resid = A %*% beta1 + B %*% wz
  dual_resid = rho * tAB %*% (wz - wz_prev)

  ## Form primal and dual tolerances.
  primal_err = ## sqrt(length(primal_resid)) * err_abs +
    err_rel * max(norm(A %*% beta1, "F"), norm(B %*% wz, "F") )
  dual_err = ## sqrt(length(dual_resid)) * err_abs +
    err_rel * norm(tA %*% U, "F")

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


##' (Helper) Check convergence for ADMM (incomplete because I couldn't figure
##' out what it should be when we are doing a matrix ADMM; will come back to
##' it).
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


##' (Helper) calculates the per-cluster objective value for the ADMM
##' @param beta p x d matrix
objective_per_cluster <- function(beta, ylist, Xa, resp, lambda, dimdat, iclust, sigma, iter,
                                  sigma_eig_by_clust=NULL){

  ## Setup
  TT = length(ylist)
  p = ncol(Xa) - 1
  mn = Xa %*% beta## matrix(b, nrow = p+1, ncol = dimdat) ## This is a (TT x dimdat)
  ## wt_resid_list = mclapply(1:TT, function(tt){
  inv.sigma = solve(sigma[iclust,,])

  ## Calculate the objective value at all times.
  terms = mclapply(1:TT, function(tt){
    y = ylist[[tt]]
    mumat = matrix(mn[tt,],
                   ncol = ncol(y),
                   nrow = nrow(y),
                   byrow = TRUE)
    wt_resid = sqrt(resp[[tt]][, iclust]) * (y - mumat)
    sums = sum(unlist(lapply(1:nrow(y), function(irow){
      t(wt_resid[irow,]) %*% inv.sigma %*% (wt_resid[irow,])
    })))

    ## Check if this is the same as the above
    browser()
    myinv_half <- sigma_eig$inverse_sigma_half
    transformed_resids = wt_resids %*% myinv_half
    sums2 = sum(transformed_resids * transformed_resids)
    stopifnot(max(abs(sums - sums2)) < 1E-8)

    return(sums)
  }, mc.cores = 8)
  all.terms = sum(unlist(terms))
  tol = 1E-8
  obj = (1/2) * sum(unlist(all.terms)) + lambda * sum(abs(beta) > tol)
  return(obj)
}


##' Calculate a specific least squares problem \min_b \|c-Db\|^2.
b_update  <- function(wvec, uw, Z, Uz, rho, yvec, D, DtD){

  cvec = c(sqrt(1/2) * yvec,
           sqrt(rho/2) * (wvec - uw/rho),
           sqrt(rho/2) * as.numeric(t(Z) - t(Uz)/rho))

  ## sol = solve(t(D) %*% D, t(D) %*% cvec)
  sol = solve(DtD, crossprod(D, cvec)) ##t(D) %*% D, t(D) %*% cvec)
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
