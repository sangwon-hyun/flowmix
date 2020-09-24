##' (Helper) Soft thresholding of |a| at radius of |b|.
##'
##' @param a Numeric vector.
##' @param b Numeric vector of same length as \code{length(b)}.
##'
##' @return Soft-thresholded version of a and b.
##'
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
##'
##' @param mat Numeric matrix.
##' @param C Radius for l2 ball projection.
##'
##' @return Projected matrix of same size as \code{mat}.
##'
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

##' (Helper) Projection into an l2 ball centered at zero, of radius C. Not used
##' now.
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
                     err_rel = 1E-4,
                     err_abs = 0
                     ){

  ## Form primal and dual residuals, and other things.
  prim1 = rbind(beta1, Xbeta1) ## (I don't like the rbind)
  prim2 = - rbind(w, Z)
  primal_resid = prim1 + prim2
  dual_resid = -rho * ((w - w_prev) + (tX %*% (Z - Z_prev)))
  tAU = Uw + tX %*% Uz

  ## Form primal and dual tolerances.
  primal_err = sqrt(length(primal_resid)) * err_abs +
    err_rel * max(norm(prim1, "F"), norm(prim2, "F"))
  dual_err = sqrt(length(dual_resid)) * err_abs +
    err_rel * norm(tAU, "F")

  ## Check convergence.
  primal_resid_size = norm(primal_resid, "F")
  dual_resid_size = norm(dual_resid, "F")
  primal_converge = ( primal_resid_size  <= primal_err )
  dual_converge = ( dual_resid_size <= dual_err )

  ## Some checks (trying to address problems with |converge|).
  assertthat::assert_that(is.numeric(primal_resid_size))
  assertthat::assert_that(is.numeric(primal_err))
  assertthat::assert_that(is.numeric(dual_resid_size))
  assertthat::assert_that(is.numeric(dual_err))

  ## return(primal_converge & dual_converge)
  converge = primal_converge & dual_converge
  if(is.na(converge)) browser()

  return(list(primal_resid = primal_resid,
              primal_err = primal_err,
              dual_resid = dual_resid,
              dual_err = dual_err,
              converge = converge))
}


##' (Helper) calculates the per-cluster objective value for the ADMM. INCOMPLETE
##' because it is super inefficient right now.
##' @param beta p x d matrix
objective_per_cluster <- function(beta, ylist, Xa, resp, lambda, N, dimdat, iclust, sigma, iter, zerothresh,
                                  first_iter, sigma_eig_by_clust=NULL, rcpp = FALSE){

  ## Setup
  TT = length(ylist)
  p = ncol(Xa) - 1
  mn = Xa %*% beta  ## VERY SLOW (600)

  ## Obtain square root of sigma.
  if(first_iter){
    sigma_half <- sigma_half_from_eig(eigendecomp_sigma_barebones(sigma[iclust,,]))
  } else {
    sigma_eig <- sigma_eig_by_clust[[iclust]]
    sigma_half <- sigma_eig$inverse_sigma_half
  }


  ## Obtain the weighted residuals once.
  ## save(ylist, TT, mn, resp, iclust, sigma_half, file="~/Desktop/rcpp-speed.Rdata")
  ## load(file="~/Desktop/rcpp-speed.Rdata", verbose=TRUE)
  if(!rcpp){
      wt_resids_list = lapply(1:TT, function(tt){
        y = ylist[[tt]]
        mumat = matrix(mn[tt,],
                       ncol = ncol(y),
                       nrow = nrow(y),
                       byrow = TRUE)
        wt_resid = sqrt(resp[[tt]][, iclust]) * (y - mumat)
      })
      transformed_resids = do.call(rbind, wt_resids_list) %*% sigma_half
      resid.prods = rowSums(transformed_resids * transformed_resids)
  }
  if(rcpp){


    ## This works and is slightly faster
    ## wt_resids_list = list()
    ## for(tt in 1:TT){
    ##   wt_resids_list[[tt]] = subtractC3(sqrt(resp[[tt]][, iclust]),
    ##                                     ylist[[tt]],
    ##                                     mn[tt,])
    ## }
    form_wt_resid <- function(tt){
      subtractC3(sqrt(resp[[tt]][, iclust]),
                 ylist[[tt]],
                 mn[tt,])
    }
    wt_resids_list = lapply(1:TT, form_wt_resid)
    transformed_resids = do.call(rbind, wt_resids_list) %*% sigma_half
    resid.prods = rowSumsC2_arma(transformed_resids)

    ## Restructuring of data.
    ## ylong = do.call(rbind, ylist) ## This can be done in \code{flowmix()}, once.
    ## ntlist = sapply(ylist, nrow)
    ## mumat = mn[rep(1:TT, ntlist), ] ## You can't get around this. ## this takes a long time.
    ## longwt = do.call(c, lapply(1:TT, function(tt){ resp[[tt]][,iclust]})) %>% sqrt()
    ## resid.prods = dothisC(longwt, ylong, mumat, sigma_half)
  }

  ## ## Calculate the inner product: resid * sigma^-1 * t(resid)
  grand.sum = sum(resid.prods)

  ## Calculate the objective value
  obj = (1/(2 * N)) * grand.sum + lambda * sum(abs(beta[-1,]) > zerothresh) ## Another bug!
  return(obj)
}


##' Calculate a specific least squares problem \deqn{\min_b \|c-Db\|^2}.
##'
##' @param wvec Numeric vector of length \code{p * dimdat}.
##' @param uw Numeric vector of length \code{p * dimdat}.
##' @param rho Positive scalar.
##' @param Z Vector of size \code{dimdat * T}.
##' @param Uz Matrix of size \code{TT} by \code{dimdat}.
##' @param yvec Vector of size \code{dimdat * T}
##' @param D Matrix of size (d(T+p) x (p+1)d)
##' @param DtDinv \deqn{(D^T D)^{-1}}
##' @param N Scalar.
##'
##' @return Least squares solution.
b_update  <- function(wvec, uw, Z, Uz, rho, yvec, D, DtDinv, N){

  cvec = c(sqrt(1/(2*N)) * yvec,
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



## Four functions to be used only within \code{admm_oneclust()}.


##' ADMM update of w vector.
##'
##' @param uw Numeric vector of length \code{p * dimdat}.
##' @param b1 Numeric vector of length \code{p * dimdat}.
##' @param rho Positive scalar.
##' @param lambda Positive scalar.
##'
##' @return Soft-thresholded version of w
wvec_update  <- function(b1, uw, lambda, rho){
  soft_thresh(b1 + uw/rho, lambda/rho)
}


##' ADMM update of Z vector.
##'
##' @param Xbeta1 Matrix of size \code{TT} by \code{dimdat}.
##' @param Uz Matrix of size \code{TT} by \code{dimdat}.
##' @param b1 Numeric vector of length \code{p * dimdat}.
##' @param C Radius, positive scalar.
##' @param rho Positive scalar.
##' @param dimdat Dimension of response data (not used).
##' @param TT Number of time points (not used).
##'
##' @return Soft-thresholded version of w
Z_update  <- function(Xbeta1, Uz, C, rho, dimdat, TT){
  mat = Xbeta1 + Uz/rho
  Z = projCmat(mat, C)
}

uw_update <- function(uw, rho, b1, wvec){
  uw + rho * (b1 - wvec)
}

Uz_update <- function(Uz, rho, Xbeta1, Z){
  Uz + rho * (Xbeta1 - Z)
}
