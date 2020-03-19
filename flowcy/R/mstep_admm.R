##' The M step of beta, using ADMM. (TODO: This should be able to use the
##' eigendecomp of the Sigmas for the objective value calculation. That is next
##' up.)
##'
##' @param maxdev The desired maximum radius of the fitted means from beta0k.
##' @param sel_coef Sparsity pattern. Only active when refit=TRUE.
##' @param sigma (numclust x dimdat x dimdat) matrix.
##' @param local_adapt TRUE if locally adaptive ADMM (LA-ADMM) is to be used. If
##'   so, \code{niter} becomes the inner number of iterations, and
##'   \code{local_adapt_niter} becomes the number of outer iterations.
##' ##' @param local_adapt_niter Number of outer iterations for LA-ADMM.
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
                            ## em_iter,

                            ## Warm startable variables
                            betas = NULL,
                            Zs = NULL,
                            wvecs=NULL,
                            uws=NULL,
                            Uzs=NULL,
                            ## End of warm startable variables

                            maxdev = NULL,
                            niter = 1E4,
                            rho = 100, ## Some default
                            err_rel = 1E-3,
                            err_abs = 0,
                            zerothresh = 1E-6,
                            plot.admm = FALSE,
                            local_adapt = FALSE,
                            local_adapt_niter = 5,
                            space = 50
                            ){

  ####################
  ## Preliminaries ###
  ####################
  TT = length(ylist)
  p = ncol(X)
  numclust = ncol(resp[[1]])
  dimdat = ncol(ylist[[1]])
  ntlist = sapply(ylist, nrow)
  N = sum(ntlist)
  Xa = cbind(1, X)
  intercept_inds = ((1:dimdat) - 1)*(p+1) + 1
  X0 = lapply(1:TT, function(tt){ diag(rep(1,dimdat)) %x% t(c(0, X[tt,,drop=TRUE]))})
  X0 = do.call(rbind, X0)
  tX = t(X)
  I_aug = make_I_aug(p, dimdat, intercept_inds)

  ##########################################
  ## Run ADMM separately on each cluster ##
  #########################################
  yhats = admm_niters = vector(length = numclust, mode = "list")
  if(first_iter) betas = vector(length = numclust, mode = "list")
  if(first_iter)  Zs =  wvecs =  uws =  Uzs = vector(length=numclust, mode="list")
  fits = matrix(NA, ncol = numclust, nrow = ceiling(niter / space))

  ## Form tilde objects for b update. Only do once!
  manip_obj = manip(ylist, Xa, resp, sigma, numclust,
                    sigma_eig_by_clust = sigma_eig_by_clust,
                    first_iter = first_iter)
  Xtildes = manip_obj$Xtildes
  yvecs = manip_obj$yvecs

  ## For every cluster, run LA-ADMM
  resid_mat_list = list()
  start.time = Sys.time()
  for(iclust in 1:numclust){
    ## cat("cluster", iclust, fill=TRUE)
    ## if(!first_iter & iclust == 2) browser()
    ## print(head(betas))

    ## Perform LA-ADMM.
    res = la_admm_oneclust(K = (if(local_adapt) local_adapt_niter else 1),
                           local_adapt = local_adapt,
                           iclust = iclust,
                           niter = niter, Xtilde = Xtildes[[iclust]], yvec = yvecs[[iclust]],
                           p = p , TT = TT, N = N, dimdat = dimdat, maxdev = maxdev, Xa = Xa,
                           rho = rho,
                           X0 = X0, I_aug = I_aug,
                           sigma = sigma,
                           intercept_inds = intercept_inds, lambda = mean_lambda,
                           resp = resp, ylist = ylist, X = X, tX = tX,
                           err_rel = err_rel,
                           err_abs = err_abs,
                           zerothresh = zerothresh,
                           sigma_eig_by_clust = sigma_eig_by_clust,
                           plot = plot.admm,
                           space = space,

                           ## Warm starts from previous *EM* iteration
                           first_iter = first_iter,
                           ## em_iter = em_iter,
                           beta = betas[[iclust]],
                           Z = Zs[[iclust]],
                           wvec=wvecs[[iclust]],
                           uw=uws[[iclust]],
                           Uz=Uzs[[iclust]])

    ## Store the results
    betas[[iclust]] = res$beta
    yhats[[iclust]] = res$yhat
    fits[,iclust] = res$fits
    admm_niters[[iclust]] = res$kk

    ## Store other things for for warmstart
    Zs[[iclust]] = res$Z
    wvecs[[iclust]] = res$wvec
    uws[[iclust]] = res$uw
    Uzs[[iclust]] = res$Uz

    resid_mat_list[[iclust]] = res$resid_mat ## temporary
  }

  ## Aggregate the yhats into one array
  yhats_array = array(NA, dim = c(TT, dimdat, numclust))
  for(iclust in 1:numclust){ yhats_array[,,iclust] = yhats[[iclust]] }

  ## Each are lists of length |numclust|.
  return(list(betas = betas,
              mns = yhats_array,
              fits = fits,
              resid_mat_list = resid_mat_list,
              admm_niters = admm_niters, ## Temporary: Seeing the number of outer
                                        ## iterations it took to converge.

              ## For warmstarts
              Zs = Zs,
              wvecs = wvecs,
              uws = uws,
              Uzs = Uzs

              ))
}


##' LA (locally adaptive) ADMM wrapper to \code{admm_oneclust()}.
la_admm_oneclust <- function(K,
                             ...){

  ## Initialize arguments for ADMM.
  args <- list(...)
  p = args$p
  TT = args$TT
  dimdat = args$dimdat
  ## em_iter = args$em_iter

    ## print(paste0("iclust=", args$iclust))
    ## print(head(args$beta[[1]]))

  ## This initialization can come from the previous *EM* iteration.
  ## print(args$first_iter)
  if(args$first_iter){
    beta = matrix(0, nrow=p+1, ncol=dimdat)
    Z = matrix(0, nrow = TT, ncol = dimdat)
    wvec = rep(0, p * dimdat)
    uw  = rep(0, p * dimdat)
    Uz = matrix(0, nrow = TT, ncol = dimdat)
    args[['beta']] <- beta
    args[['Z']] <- Z
    args[['wvec']] <- wvec
    args[['uw']] <- uw
    args[['Uz']] <- Uz
  }

  fits = c()
  cols = c()
  objectives = c()

  outer_converge <- function(objectives){
    ## print(objectives)
    consec = 4
    if(length(objectives) < consec){
      return(FALSE)
    } else {
      mytail = tail(objectives, consec)
      rel_diffs = mytail[1:(consec-1)]/mytail[2:consec]
      ## print(rel_diffs)
      return(all(abs(rel_diffs) - 1 < 1E-3))
    }
  }

  ## Run ADMM repeatedly with (1) double rho, and (2) previous b
  for(kk in 1:K){
    ## printprogress(kk, K, "outer admm", fill=TRUE)
    if(kk>1){
    args[['beta']] <- beta
    args[['Z']] <- Z
    args[['wvec']] <- wvec
    args[['uw']] <- uw
    args[['Uz']] <- Uz
    }

    ## Run ADMM
    if(kk > 1) args[['rho']] <- rho
    ## args[['beta']] <- beta
    ## args[['Z']] <- Z
    ## args[['wvec']] <- wvec
    ## args[['uw']] <- uw
    ## args[['Uz']] <- Uz


    ## 1. Old call of function, which is slow and memory intensive.
    ## res = do.call(admm_oneclust, args)

    ## 2. New way no.1: Courtesy of Hadley Wickham
    argn <- lapply(names(args), as.name)
    names(argn) <- names(args)
    call <- as.call(c(list(as.name("admm_oneclust")), argn))
    res = eval(call, args)

    ## 3. New way no.2: Courtesy of rmisc.
    ## res = rlang::invoke("admm_oneclust", args)

    ## See if outer iterations should terminate
    objectives = c(objectives, res$fit)
      ## if(em_iter>2) plot(fits, type='o', col=cols)
    if(outer_converge(objectives)){
      ## cat("Outer ADMM converged with K=", kk, fill=TRUE)

      ## ylimlist = list(c(3.192, 3.291), c(2.918, 3.027), c(5.27,5.92), c(0.835,0.87))
      ## ylimlist = lapply(ylimlist, function(a) a*c(0.8, 1.2))
      ## xlimlist = list(c(0,120), c(0,100), c(0,120), c(0, 120))
      ## if(em_iter > 2) plot(fits, type='o', col=cols)
                         ## xlim=xlimlist[[args$iclust]],
                         ## ylim=ylimlist[[args$iclust]],
                         ## main = "Yes EM warmstart")##, main = paste0("rho=", signif(args$rho,3)))
      break
    }
    if(res$converge){
      ## cat("converged with K=", kk, fill=TRUE)
      break
    }

    ## Update some parameters
    rho = args[['rho']]
    rho = rho * 2
    beta = res$beta
    Z = res$Z
    wvec = res$wvec
    uw = res$uw
    Uz = res$Uz

  }
  ## Sys.sleep(5) ## Temporary
  ## if(!res$converge) print("Didn't converge at all") ## Change to warning

  ## Record how long the admm took; in terms of # iterations.
  res$kk = kk

  return(res)
}


##' ADMM of one cluster.
##'
##' @param local_adapt TRUE if locally adaptive ADMM is to be used.
##'
##' @return List containing |beta|, |yhat|, |resid_mat|, |fits|.
admm_oneclust <- function(iclust, niter, Xtilde, yvec, p,
                          TT, N, dimdat, maxdev, Xa, rho,
                          X0, I_aug,
                          intercept_inds, lambda,
                          resp, ylist, X, tX, err_rel, err_abs,
                          zerothresh,
                          ## Warm startable variables
                          beta,
                          Z,
                          wvec,
                          uw,
                          Uz,
                          first_iter,## Not used
                          ## em_iter,
                          ## End of warm startable variables
                          local_adapt,
                          sigma_eig_by_clust,
                          sigma,
                          space = 20,
                          plot = FALSE){

  resid_mat = matrix(NA, nrow = ceiling(niter/5), ncol = 4)
  colnames(resid_mat) = c("primresid", "primerr", "dualresid", "dualerr")

  ## Prepare an object for b_update()
  Dfirst = sqrt(1/(2*N)) * Xtilde
  Drest = rbind(sqrt(rho/2) * I_aug,
                sqrt(rho/2) * X0)
  D = rbind(Dfirst, Drest)
  DtD = crossprod(Dfirst, Dfirst) + crossprod(Drest, Drest)
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
  ## Z = matrix(0, nrow = TT, ncol = dimdat)
  ## wvec = uw = rep(0, p * dimdat)
  ## Uz = matrix(0, nrow = TT, ncol = dimdat)
  C = maxdev
  fits = rep(NA, ceiling(niter/space))
  converge = FALSE

  for(iter in 1:niter){
    ## printprogress(iter, niter, "inner admm")

    b <- b_update(wvec, uw, Z, Uz, rho, yvec, D, DtDinv, N)
    b1 = b[-intercept_inds]
    b0 = b[intercept_inds]
    beta = matrix(b, nrow = p+1)
    beta1 = matrix(b1, nrow = p)
    Xbeta1 = X %*% beta1

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

    ## 3. Check convergence
    if( iter > 1  & iter %% 5 == 0){## & !local_adapt){

      ## Calculate convergence criterion
      obj = converge(beta1, rho, w, Z, w_prev, Z_prev, Uw, Uz, tX = tX,
                     Xbeta1 = Xbeta1, err_rel = err_rel,
                     err_abs = err_abs)
      jj = (iter/ 5)
      resid_mat[jj,] = c(norm(obj$primal_resid, "F"),
                         obj$primal_err,
                         norm(obj$dual_resid,"F"),
                         obj$dual_err)

      if(obj$converge){
        print(paste('Converged! in', iter, 'out of ', niter, 'steps!')) ## ADMM debug code
        converge = TRUE
        break
      }
    }


    ## 3. Calculate objective values for this cluster.
    w_prev = w
    Z_prev = Z

    ## ## 4. Calculate things related to convergence (Slow).
    if(iter %% space == 0 ){
      ii = iter / space
      ## if(em_iter>2){
      ## fits[ii] = objective_per_cluster(beta, ylist, Xa, resp, lambda,
      ##                                  N, dimdat, iclust, sigma, iter,
      ##                                  zerothresh,
      ##                                  TRUE) ## Just flagging first_iter=TRUE for now.
      ## }
    }
  }

  ## browser()
        ## plot(fits, type = 'l', main = paste("cluster", iclust))

  ## if(!obj$converge){
  ##       plot((resid_mat[,"primresid"]), type = 'l', main = paste("Primal resid, cluster", iclust),
  ##            ylim = range((resid_mat[,c("primresid", "primerr")]), na.rm=TRUE))
  ##       lines((resid_mat[,"primerr"]), type = 'l', col='red', lwd=2)

  ##       plot(log(resid_mat[,"dualresid"]), type = 'l', main = paste("Dual resid, cluster", iclust),
  ##            ylim = range(log(resid_mat[,c("dualresid", "dualerr")]), na.rm=TRUE))
  ##       lines(log(resid_mat[,"dualerr"]), type = 'l', col='red', lwd=2)
  ## }

  ## Temporary print message
  ## if(!(obj$converge)){
  ##   print(paste("Didn't converge in", niter, 'steps!'))
  ## }

  ## Gather results.
  beta[-1,] = w
  beta[which(abs(beta) < zerothresh, arr.ind = TRUE)] = 0
  yhat = Xa %*% beta

  fit = objective_per_cluster(beta, ylist, Xa, resp, lambda, N, dimdat,
                              iclust, sigma, iter, zerothresh,
                              is.null(sigma_eig_by_clust), sigma_eig_by_clust)

  return(list(beta = beta, yhat = yhat, resid_mat = resid_mat, fits = fits, converge=converge,
              fit=fit,
              ## Other variables to return.
              Z = Z,
              wvec = wvec,
              uw = uw,
              Uz = Uz))
}
