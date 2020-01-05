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
                            rho = 0.1, ## Some default
                            err_rel = 1E-3,
                            converge_fast=TRUE ## temporary
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

  if(converge_fast){
    tX = t(X)
  } else {
    ## Prepare a few objects for converge()
    p = ncol(X)
    TT = nrow(X)
    A = rbind(diag(rep(1,p)),
              X)
    tA = t(A)
    B = Matrix::bdiag(-diag(rep(1,p)),
                      -diag(rep(1,TT)))
    B = as.matrix(B)
    tAB = tA %*% B
  }

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

  ## Prepare a few objects for converge()
  p = ncol(X)
  TT = nrow(X)

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

    ## Prepare an object for b_update()
    D = rbind(sqrt(1/2) * Xtildes[[iclust]],
              sqrt(rho/2) * I_aug,
              sqrt(rho/2) * X0) ## This can be moved out
    DtD = crossprod(D, D)

    for(iter in 1:niter){

      b <- b_update(wvec, uw, Z, Uz, rho, yvecs[[iclust]], D, DtD)
      b1 = b[-intercept_inds]
      b0 = b[intercept_inds]
      beta = matrix(b, nrow = p+1)
      beta1 = matrix(b1, nrow = p)
      Xbeta1 = X %*% beta1

      Z_update  <- function(Xbeta1, Uz, C, rho, dimdat, TT){
        mat = Xbeta1 + Uz/rho
        Z = projCmat(mat, C)
      }
      Z <- Z_update(Xbeta1, Uz, C, rho, dimdat, TT)

      wvec_update  <- function(b1, uw, lambda, rho){
        soft_thresh(b1 + uw/rho, lambda/rho)
      }
      wvec <- wvec_update(b1, uw, lambda, rho)
      w <- matrix(wvec, nrow = p, byrow=FALSE)

      Uz_update <- function(Uz, rho, Xbeta1, Z){
        Uz + rho * (Xbeta1 - Z)
      }
      Uz <- Uz_update(Uz, rho, Xbeta1, Z)

      uw_update <- function(uw, rho, b1, wvec){
        uw + rho * (b1 - wvec)
      }
      uw <- uw_update(uw, rho, b1, wvec)
      Uw <- matrix(uw, nrow = p, byrow=FALSE)

      ## 3. Check convergence
      if( iter > 1  & iter %% 5 == 0 ){
        if(converge_fast){
          obj = converge(beta1, rho, w, Z, w_prev, Z_prev, Uw, Uz, tX = tX, Xbeta1 = Xbeta1, err_rel = err_rel)

        } else{
          obj = converge_old(beta1, X, rho, w, Z, w_prev, Z_prev, Uw, Uz, A=A, B=B, tA=tA, tAB = tAB, err_rel = err_rel)
        }
        resid_mat[iter-1,] = c(norm(obj$primal_resid, "F"), obj$primal_err, norm(obj$dual_resid,"F"), obj$dual_err)
        if(obj$converge){  break }
      }
      w_prev = w
      Z_prev = Z

      ## 3. Calculate objective values for this cluster.
      ## if(iter %% 20 == 0){
      ##   fits[iter, iclust] = objective_per_cluster(beta, ylist, Xa, resp, lambda,
      ##                                                                dimdat, iclust, sigma, iter)
      ## }
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
