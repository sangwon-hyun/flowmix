##' The M step of beta, using ADMM. (TODO: This should be able to use the
##' eigendecomp of the Sigmas for the objective value calculation. That is next
##' up.)
##'
##' @param maxdev The desired maximum radius of the fitted means from beta0k.
##' @param sel_coef Sparsity pattern. Only active when refit=TRUE.
##' @param sigma (numclust x dimdat x dimdat) matrix.
##' @param local_adapt TRUE if locally adaptive ADMM (LA-ADMM) is to be used. If
##'   so, \code{niter} becomes the inner number of iterations, and
##'   \code{local_adapt_niter} becomes the number of outer iterations. Also, if
##'   TRUE, \code{err_rel} is not used.
##' @param local_adapt_niter Number of outer iterations for LA-ADMM.
##'
##' @return Result of M step; a |numclust| length list of (p+1)x(d) matrices,
##'   each containing the estimated coefficients for the mean estimation.
Mstep_beta_admm <- function(resp,
                            ylist,
                            X,
                            mean_lambda = 0,
                            sigma,
                            sigma_eig_by_clust = NULL,
                            first_iter = TRUE,
                            maxdev = NULL,
                            niter = 1E4,
                            rho = 100, ## Some default
                            err_rel = 1E-3,
                            zerothresh = 1E-6,
                            plot.admm = FALSE,
                            local_adapt = FALSE,
                            local_adapt_niter = 5
                            ){

  ####################
  ## Preliminaries ###
  ####################
  TT = length(ylist)
  p = ncol(X)
  numclust = ncol(resp[[1]])
  dimdat = ncol(ylist[[1]])
  ntlist = sapply(ylist, nrow)
  Xa = cbind(1, X)
  intercept_inds = ((1:dimdat) - 1)*(p+1) + 1
  X0 = lapply(1:TT, function(tt){ diag(rep(1,dimdat)) %x% t(c(0, X[tt,,drop=TRUE]))})
  X0 = do.call(rbind, X0)
  lambda = mean_lambda
  tX = t(X)
  I_aug = make_I_aug(p, dimdat, intercept_inds)

  ##########################################
  ## Run ADMM separately on each cluster ##
  #########################################
  betas = yhats = vector(length = numclust, mode = "list")
  fits = matrix(NA, ncol = numclust, nrow = ceiling(niter / 20))

  ## 1. Form tilde objects for b update. Only do once!
  manip_obj = manip(ylist, Xa, resp, sigma, numclust,
                    sigma_eig_by_clust = sigma_eig_by_clust,
                    first_iter = first_iter)
  Xtildes = manip_obj$Xtildes
  yvecs = manip_obj$yvecs

  ## Temporary
  resid_mat_list = list()
  start.time = Sys.time()
  for(iclust in 1:numclust){

    ## Perform ADMM once.
    res = la_admm_oneclust(K = (if(local_adapt) local_adapt_niter else 1),
                           local_adapt = local_adapt,
                           iclust = iclust,
                           niter = niter, Xtilde = Xtildes[[iclust]], yvec = yvecs[[iclust]],
                           p =p , TT = TT, dimdat = dimdat, maxdev = maxdev, Xa = Xa,
                           rho = rho,
                           X0, I_aug,
                           intercept_inds = intercept_inds, lambda = lambda,
                           resp = resp, ylist = ylist, tX = tX,
                           err_rel = err_rel,
                           zerothresh = zerothresh,
                           plot = plot.admm)

    ## Store the results (only b)
    betas[[iclust]] = res$beta
    yhats[[iclust]] = res$yhat
    fits[,iclust] = res$fits
    resid_mat_list[[iclust]] = res$resid_mat ## temporary
  }

  ## Aggregate the yhats into one array
  yhats_array = array(NA, dim = c(TT, dimdat, numclust))
  for(iclust in 1:numclust){ yhats_array[,,iclust] = yhats[[iclust]] }

  ## Each are lists of length |numclust|.
  return(list(beta = betas,
              mns = yhats_array,
              fits = fits,
              resid_mat_list = resid_mat_list
              ))
}


##' LA (locally adaptive) ADMM wrapper to \code{admm_oneclust()}.
la_admm_oneclust <- function(K = 5,
                             ...){

  ## Initialize arguments for ADMM.
  args <- list(...)
  p = args$p
  dimdat = args$dimdat
  beta = matrix(0, nrow=p+1, ncol=dimdat)
  args[['beta']] <- beta

  ## Run ADMM repeatedly with (1) double rho, and (2) previous b
  for(kk in 1:K){
    ## printprogress(kk, K, "outer admm", fill=TRUE)

    ## Run ADMM
    if(kk > 1) args[['rho']] <- rho
    args[['beta']] <- beta
    res = do.call(admm_oneclust, args)

    ## Update some parameters
    rho = args[['rho']]
    rho = rho * 2
    beta = res$beta
  }

  return(res)
}


##' ADMM of one cluster.
##'
##' @param local_adapt TRUE if locally adaptive ADMM is to be used.
##'
##' @return List containing |beta|, |yhat|, |resid_mat|, |fits|.
admm_oneclust <- function(iclust, niter, Xtilde, yvec, p,
                          TT, dimdat, maxdev, Xa, rho,
                          X0, I_aug,
                          intercept_inds, lambda,
                          resp, ylist, tX, err_rel,
                          zerothresh,
                          beta,
                          local_adapt,
                          plot = FALSE){

  resid_mat = matrix(NA, nrow = ceiling(niter/5), ncol = 4)
  colnames(resid_mat) = c("primresid", "primerr", "dualresid", "dualerr")

  ## Prepare an object for b_update()
  Drest = rbind(sqrt(rho/(2*TT)) * I_aug,
                sqrt(rho/(2*TT)) * X0)
  Drest_square = crossprod(Drest, Drest)
  Dfirst = sqrt(1/2) * Xtilde
  D = rbind(Dfirst, Drest)
  DtD = crossprod(Dfirst, Dfirst) + Drest_square

  ## DtD = crossprod(D, D)
  DtDinv = chol2inv(chol(DtD))

  ###############################
  ## Initialize the variables ###
  ###############################
  b = as.numeric(beta)
  b1 = b[-intercept_inds]
  b0 = b[intercept_inds]
  beta1 = matrix(b1, nrow = p)
  Xbeta1 = X %*% beta1
  Z = matrix(0, nrow = TT, ncol = dimdat)
  wvec = uw = rep(0, p * dimdat)
  Uz = matrix(0, nrow = TT, ncol = dimdat)
  C = maxdev
  fits = rep(NA, ceiling(niter/20))

  for(iter in 1:niter){
    ## printprogress(iter, niter, "inner admm")

    ## If locally adaptive ADMM is used, first iteration is skipped because a
    ## warmstart beta has been provided.
    if( iter > 1){
      b <- b_update(wvec, uw, Z, Uz, rho, yvec, D, DtDinv)
      b1 = b[-intercept_inds]
      b0 = b[intercept_inds]
      beta = matrix(b, nrow = p+1)
      beta1 = matrix(b1, nrow = p)
      Xbeta1 = X %*% beta1
    }

    wvec_update  <- function(b1, uw, lambda, rho){
      soft_thresh(b1 + uw/rho, lambda/rho)
    }
    wvec <- wvec_update(b1, uw, lambda, rho)
    w <- matrix(wvec, nrow = p, byrow=FALSE)


    Z_update  <- function(Xbeta1, Uz, C, rho, dimdat, TT){
      mat = Xbeta1 + Uz/rho
      Z = projCmat(mat, C)
    }
    Z <- Z_update(Xbeta1, Uz, C, rho, dimdat, TT)

    uw_update <- function(uw, rho, b1, wvec){
      uw + rho * (b1 - wvec)
    }
    uw <- uw_update(uw, rho, b1, wvec)
    Uw <- matrix(uw, nrow = p, byrow=FALSE)

    Uz_update <- function(Uz, rho, Xbeta1, Z){
      Uz + rho * (Xbeta1 - Z)
    }
    Uz <- Uz_update(Uz, rho, Xbeta1, Z)
    ## print(max(as.numeric(Uz)))

    ## 3. Check convergence
    if( iter > 1  & iter %% 5 == 0 & !local_adapt){
      obj = converge(beta1, rho, w, Z, w_prev, Z_prev, Uw, Uz, tX = tX,
                     Xbeta1 = Xbeta1, err_rel = err_rel)
      ## ii = iter-1
      jj = (iter/ 5)
      resid_mat[jj,] = c(norm(obj$primal_resid, "F"),
                         obj$primal_err,
                         norm(obj$dual_resid,"F"),
                         obj$dual_err)
      ## Temporary print message
      if(obj$converge){
        print(paste('converged! in', iter, 'out of ', niter, 'steps!'))
        break
      }
    }


    ## 3. Calculate objective values for this cluster.
    w_prev = w
    Z_prev = Z

    ## 4. Calculate things related to convergence.
    ## if(iter %% 20 == 0 ){
    ##   ii = iter / 20
    ##   fits[ii] = objective_per_cluster(beta, ylist, Xa, resp, lambda,
    ##                                    dimdat, iclust, sigma, iter,
    ##                                    zerothresh,
    ##                                    FALSE,
    ##                                    sigma_eig_by_clust)
    ##   ## if(plot){
    ##   ##   par(mfrow=c(1,3))
    ##   ##   plot(fits, type = 'l', main = paste("cluster", iclust))
    ##   ##   plot((resid_mat[,"primresid"]), type = 'l', main = paste("Primal resid, cluster", iclust),
    ##   ##        ylim = range((resid_mat[,c("primresid", "primerr")]), na.rm=TRUE))
    ##   ##   lines((resid_mat[,"primerr"]), type = 'l', col='red', lwd=2)
    ##   ##   plot((resid_mat[,"dualresid"]), type = 'l', main = paste("Dual resid, cluster", iclust),
    ##   ##        ylim = range((resid_mat[,c("dualresid", "dualerr")]), na.rm=TRUE))
    ##   ##   lines((resid_mat[,"dualerr"]), type = 'l', col='red', lwd=2)
    ##   ## }
    ## }
  }

  ## Temporary print message
  if(!local_adapt ){
    if(!(obj$converge)){
      print(paste("Didn't converge in", niter, 'steps!'))
    }
  }

  ## Gather results.
  beta[-1,] = w
  beta[which(abs(beta) < zerothresh, arr.ind = TRUE)] = 0
  yhat = Xa %*% beta

  return(list(beta = beta, yhat = yhat, resid_mat = resid_mat, fits = fits))
}
