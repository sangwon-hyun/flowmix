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
                            space = 50,
                            rcpp = FALSE
                            ){

  ####################
  ## Preliminaries ###
  ####################
  TT = length(ylist)
  p = ncol(X)
  numclust = ncol(resp[[1]])
  dimdat = ncol(ylist[[1]])
  ntlist = sapply(ylist, nrow)
  ## N = sum(ntlist) ## OLD BUG!!
  resp.sum = t(sapply(resp, colSums)) ## (T x numclust)
  resp.sum = as.matrix(resp.sum)
  N = sum(resp.sum) ## NEW (make more efficient, later)

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
  if(first_iter)  Zs =  wvecs =  uws =  Uzs = vector(length = numclust, mode = "list")
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

  ## ## Begin temporary
  ## Prepare an object for b_update()
  Dfirst = sqrt(1/(2*N)) * Xtildes[[iclust]]
  Drest = rbind(sqrt(rho/2) * I_aug,
                sqrt(rho/2) * X0)
  D = rbind(Dfirst, Drest)
  DtD = crossprod(Dfirst, Dfirst) + crossprod(Drest, Drest)
  DtDinv = chol2inv(chol(DtD))
  Dobj = DtDinv %*% t(D)
  ## ## End temporary


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
                           Dobj = Dobj,## Temporary
                           sigma_eig_by_clust = sigma_eig_by_clust,
                           plot = plot.admm,
                           space = space,
                           rcpp = rcpp,

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
    ## fits[,iclust] = res$fits
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
              admm_niters = admm_niters, ## Temporary: Seeing the number of
                                         ## outer iterations it took to
                                         ## converge.

              ## For warmstarts
              Zs = Zs,
              wvecs = wvecs,
              uws = uws,
              Uzs = Uzs

              ))
}


##' LA (locally adaptive) ADMM wrapper to \code{admm_oneclust()}.
##'
##' @param K Number of outer iterations.
##'
la_admm_oneclust <- function(K,
                             ...){

  ## Initialize arguments for ADMM.
  args <- list(...)
  p = args$p
  TT = args$TT
  dimdat = args$dimdat

  ## This initialization can come from the previous *EM* iteration.
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

  cols = c()
  objectives = c()

  ## Temporary
  all_fits = c()
  all_cols = c()

  ## Run ADMM repeatedly with (1) double rho, and (2) previous b
  for(kk in 1:K){
    if(kk > 1){
      args[['beta']] <- beta
      args[['Z']] <- Z
      args[['wvec']] <- wvec
      args[['uw']] <- uw
      args[['Uz']] <- Uz
    }

    ## Run ADMM
    if(kk > 1) args[['rho']] <- rho

    ## Call main function
    argn <- lapply(names(args), as.name)
    names(argn) <- names(args)
    call <- as.call(c(list(as.name("admm_oneclust")), argn))
    res = eval(call, args)


    ## ## Temporary (uncomment for plotting objectives)
    ## fits = as.numeric(na.omit(res$fits))
    ## all_fits = c(all_fits, fits)
    ## all_cols = c(all_cols, rep(kk, length(fits)))
    ## print(all_fits)
    ## if(length(all_fits)!=0){
    ##   plot(all_fits %>% log10(), col=all_cols, type='o', lwd=2, main = paste0("outer iter = ", kk))
    ## }
    ## ## End temporary


    ## See if outer iterations should terminate
    objectives = c(objectives, res$fit)
    if(outer_converge(objectives) | res$converge){
      break
    }

    ## Update some parameters; double the rho value
    rho = args$rho * 2
    beta = res$beta
    Z = res$Z
    wvec = res$wvec
    uw = res$uw
    Uz = res$Uz

  }
  ## if(!res$converge) warning("Didn't converge at all")

  ## Record how long the admm took; in terms of # iterations.
  res$kk = kk

  return(res)
}

##' Check if 4 consecutive objective values are sufficiently close to 1.
##'
##' @param objectives Numeric vector.
##'
##' @return Boolean.
outer_converge <- function(objectives){
  consec = 4
  if(length(objectives) < consec){
    return(FALSE)
  } else {
    mytail = utils::tail(objectives, consec)
    rel_diffs = mytail[1:(consec-1)]/mytail[2:consec]
    return(all(abs(rel_diffs) - 1 < 1E-3))
  }
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
                          Dobj, ## Temporary
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
                          plot = FALSE,
                          rcpp = FALSE){

  resid_mat = matrix(NA, nrow = ceiling(niter/5), ncol = 4)
  colnames(resid_mat) = c("primresid", "primerr", "dualresid", "dualerr")

  ## Prepare an object for b_update()
  Dfirst = sqrt(1/(2*N)) * Xtilde
  Drest = rbind(sqrt(rho/2) * I_aug,
                sqrt(rho/2) * X0)
  D = rbind(Dfirst, Drest)
  DtD = crossprod(Dfirst, Dfirst) + crossprod(Drest, Drest)
  DtDinv = chol2inv(chol(DtD))
  Dobj = DtDinv %*% t(D)

  ###############################
  ## Initialize the variables ###
  ###############################
  b = as.numeric(beta)
  b1 = b[-intercept_inds]
  b0 = b[intercept_inds]
  beta1 = matrix(b1, nrow = p)
  Xbeta1 = X %*% beta1
  C = maxdev
  fits = rep(NA, ceiling(niter/space))
  converge = FALSE

  for(iter in 1:niter){

    if(!rcpp){
      b <- b_update(wvec, uw, Z, Uz, rho, yvec, D, DtDinv, N)
    }
    if(rcpp){
      cvec3_el = as.numeric(t(Z - Uz/rho)) ## This is still slow; make it faster!
      b <- b_updateC(wvec, uw, rho, cvec3_el, yvec, D, DtDinv, N)
    }
    b1 = b[-intercept_inds]
    b0 = b[intercept_inds]
    beta = matrix(b, nrow = p+1)
    beta1 = matrix(b1, nrow = p)
    Xbeta1 = X %*% beta1

    ## Update 1
    if(!rcpp)  wvec <- wvec_update(b1, uw, lambda, rho)
    if(rcpp)   wvec <- wvec_updateC(b1, uw, lambda, rho)
    w <- matrix(wvec, nrow = p, byrow=FALSE)

    ## Update 2
    if(!rcpp) Z <- Z_update(Xbeta1, Uz, C, rho)##, dimdat, TT)
    if(rcpp)  Z <- Z_updateC(Xbeta1, Uz, C, rho, dimdat, TT)

    ## Update 3
    uw <- uw_update(uw, rho, b1, wvec)
    Uw <- matrix(uw, nrow = p, byrow=FALSE)

    ## Update 4
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
        converge = TRUE
        break
      }
    }


    ## 3. Calculate objective values for this cluster.
    w_prev = w
    Z_prev = Z

    ## ## Temporary (uncomment for plotting objectives)
    ## ## 4. Calculate things related to convergence (Slow).
    ## if(iter %% space == 0 ){
    ##   ii = iter / space
    ##   fits[ii] = objective_per_cluster(beta, ylist, Xa, resp, lambda, N, dimdat,
    ##                               iclust, sigma, iter, zerothresh,
    ##                               is.null(sigma_eig_by_clust), sigma_eig_by_clust)
    ## }
    ## ## End of temporary
  }

  ## plot(fits, type='l', title=iclust)
  ## Sys.sleep(10)

  ## Gather results.
  beta[-1,] = w
  beta[which(abs(beta) < zerothresh, arr.ind = TRUE)] = 0
  yhat = Xa %*% beta

  fit = objective_per_cluster(beta, ylist, Xa, resp, lambda, N, dimdat,
                              iclust, sigma, iter, zerothresh,
                              is.null(sigma_eig_by_clust), sigma_eig_by_clust,
                              rcpp = rcpp)

  return(list(beta = beta,
              yhat = yhat,
              resid_mat = resid_mat,
              fits = fits,
              converge = converge,
              fit = fit,
              ## Other variables to return.
              Z = Z,
              wvec = wvec,
              uw = uw,
              Uz = Uz))
}
