##' The M step of beta, using ADMM. (TODO: This should be able to use the
##' eigendecomp of the Sigmas for the objective value calculation. That is next
##' up.)
##'
##' @param maxdev The desired maximum radius of the fitted means from beta0k.
##' @param sel_coef Sparsity pattern. Only active when refit=TRUE.
##' @param sigma (numclust x dimdat x dimdat) matrix.
##'
##' @return Result of M step; a |numclust| length list of (p+1)x(d) matrices,
##'   each containing the estimated coefficients for the mean estimation.
Mstep_beta_admm <- function(resp, ylist, X, mean_lambda = 0, sigma, numclust,
                            maxdev = NULL,
                            niter = 100,
                            rho = 0.1 ## Some default
                            ){

  ####################
  ## Preliminaries ###
  ####################
  TT = length(ylist)
  numclust = ncol(resp[[1]])
  dimdat = ncol(ylist[[1]])
  ntlist = sapply(ylist, nrow)
  p = ncol(X)
  Xa = cbind(1, X)
  intercept_inds = ((1:dimdat) - 1)*(p+1) + 1
  X0 = lapply(1:TT, function(tt){ diag(rep(1,dimdat)) %x% t(c(0, X[tt,,drop=TRUE]))})
  X0 = do.call(rbind, X0)
  lambda = mean_lambda
  I_aug = make_I_aug(p, dimdat, intercept_inds)

  ##########################################
  ## Run ADMM separately on each cluster ##
  #########################################
  betas = yhats = vector(length = numclust, mode = "list")
  fits = matrix(0, ncol = numclust, nrow = niter)

  ###############################
  ## Initialize the variables ###
  ###############################
  b = rep(0, dimdat*(p+1))
  Z = matrix(0, nrow = TT, ncol = dimdat)
  wvec = uw = rep(0, p * dimdat)
  Uz = matrix(0, nrow = TT, ncol = dimdat)
  C = maxdev

  ## 1. Form tilde objects for b update. Only do once!
  manip_obj = manip(ylist, Xa, resp, sigma, numclust,
                    first_iter = TRUE) ## todo: consider updating this.
  Xtildes = manip_obj$Xtildes
  yvecs = manip_obj$yvecs


  ## Temporary
  resid_mat_list = list()
  start.time = Sys.time()
  for(iclust in 1:numclust){
    resid_mat = matrix(0, nrow = niter, ncol = 4)
    ## printprogress(iclust, numclust, "cluster", start.time = start.time, fill=TRUE)

    for(iter in 1:niter){
      ## printprogress(iter, niter, "iteration")

      b = b_update(wvec, uw, Z, Uz, X0, rho, Xtildes[[iclust]], yvecs[[iclust]], I_aug)
      b1 = b[-intercept_inds]
      b0 = b[intercept_inds]
      beta = matrix(b, nrow = p+1)
      beta1 = matrix(b1, nrow = p)
      assert_that(all.equal(as.numeric(I_aug %*% b ), b1) == TRUE) ## Check once

      Z_update  <- function(beta1, X, Uz, C, rho){
        mat = X %*% beta1 + Uz/rho
        t(apply(mat, 1, projC, C)) ## TODO: improve.
      }
      Z <- Z_update(beta1, X, Uz, C, rho)

      wvec_update  <- function(b1, uw, lambda, rho){
        soft_thresh(b1 + uw/rho, lambda/rho)
      }
      wvec <- wvec_update(b1, uw, lambda, rho)
      w <- matrix(wvec, nrow = p, byrow=FALSE)

      Uz_update <- function(Uz, rho, beta1, X, Z){
        Uz + rho * (X %*% beta1 - Z)
      }
      Uz <- Uz_update(Uz, rho, beta1, X, Z)

      uw_update <- function(uw, rho, b1, wvec){
        uw + rho * (b1 - wvec)
      }
      uw <- uw_update(uw, rho, b1, wvec)
      Uw <- matrix(uw, nrow = p, byrow=FALSE)

      ## 3. Check convergence
      if( iter > 1 ){
        ## if( converge(beta1, X, rho, w, Z, w_prev, Z_prev, Uw, Uz) ) next
        obj = converge(beta1, X, rho, w, Z, w_prev, Z_prev, Uw, Uz)
        resid_mat[iter-1,] = c(norm(obj$primal_resid, "F"), obj$primal_err, norm(obj$dual_resid,"F"), obj$dual_err)
        if(obj$converge) break
      }
      w_prev = w
      Z_prev = Z

      ## 3. Calculate objective values for this cluster.
      if(iter %% 20 == 0){
        fits[iter, iclust] = objective_per_cluster(beta, ylist, Xa, resp, lambda,
                                                                     dimdat, iclust, sigma, iter)
      }
    }

    ## Store the results (only b)
    yhat = Xa %*% beta
    betas[[iclust]] = beta
    yhats[[iclust]] = yhat
    resid_mat_list[[iclust]] = resid_mat ## temporary
  }

  ## Aggregate the yhats into one array
  yhats_array = array(NA, dim = c(TT, dimdat, numclust))
  for(iclust in 1:numclust){ yhats_array[,,iclust] = yhats[[iclust]] }

  ## Each are lists of length |numclust|.
  return(list(beta = betas,
              mns = yhats_array,
              fits = fits, ## Temporary
              resid_mat_list = resid_mat_list ## Temporary
              ))
}
