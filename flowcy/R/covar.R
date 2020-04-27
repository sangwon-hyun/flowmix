##' Main function for covariate EM. Does covariate EM with |nrep| restarts (5 by
##' default).
##' @param ... Arguments for \code{covarem_once()}.
##' @param nrep Number of restarts.
##' @return The |covarem| class object that had the best likelihood.
##' @export
covarem <- function(..., nrep = 5){

  ## Don't do many restarts if warmstart-able mean exists
  dots <- list(...)
  if(!is.null(dots$mn)) nrep = 1

  ## Do |nrep| restarts.
  reslist = list()
  for(itrial in 1:nrep){
    reslist[[itrial]] = covarem_once(...)
  }

  ## Pick the best one and return
  objlist = lapply(reslist, function(res){ res$obj[-1]})
  ii = which.min(sapply(objlist, min))
  return(reslist[[ii]])
}

##' Main function for covariate EM.
##'
##' @param ylist T-length list each containing response matrices of size (nt x
##'   3), which contains coordinates of the 3-variate particles, organized over
##'   time (T) and with (nt) particles at every time.
##' @param X Matrix of size (T x p+1)
##' @param pie.list (T by K)
##' @param mean_lambda lambda for lasso for the mean.
##' @param pie_lambda lambda for lasso for pie.
##' @param refit (experimental), defaults to FALSE. If TRUE, then the refitted
##'   non-regularized solutions (with only a user-specified set of active
##'   coefficients, coded in the argument \code{sel_coef}) are calculated.
##' @param sel_coef (Experimental feature) Boolean matrices of the same
##'   structure as beta and alpha, whose TRUE entries are the active
##'   coefficients to be refitted in a nonregularized way.
##' @param tol_em Relative tolerance for EM convergence. Defaults to 1E-4.
##' @param zero_stabilize Defaults to TRUE, in which case the EM is only run
##'   until the zero pattern in the coefficients stabilize.
##'
##' @return List containing fitted parameters and means and mixture weights,
##'   across algorithm iterations. beta is a (p+1 x 3 x dimdat) array. Alpha is
##'   a (dimdat x (p+1)) array.
##'
##' @export
covarem_once <- function(ylist, X,
                         countslist = NULL,
                         numclust, niter = 1000,
                         mn = NULL, pie_lambda,
                         mean_lambda, verbose = FALSE,
                         sigma.fac = 1, tol_em = 1E-4,
                         refit = FALSE, ## EXPERIMENTAL FEATURE.
                         sel_coef = NULL,
                         maxdev = NULL,
                         manual.bin = FALSE,
                         manual.grid = NULL,
                         countslist_overwrite = NULL,
                         zero_stabilize  = FALSE,
                         ## Temporary
                         ridge = FALSE,
                         ridge_lambda = 0,
                         ## End of temporary
                         plot = FALSE,
                         plotdir = "~/Desktop",
                         init_mn_flatten = FALSE,
                         ## beta Mstep (CVXR) settings
                         mstep_cvxr_ecos_thresh = 1E-8,
                         mstep_cvxr_scs_eps = 1E-5,
                         zerothresh = 1E-6,
                         ## beta Mstep (ADMM) settings
                         admm = TRUE,
                         admm_rho = 0.01,
                         admm_err_rel = 1E-3,
                         admm_err_abs = 1E-4,
                         ## beta M step (Locally Adaptive ADMM) settings
                         admm_local_adapt = TRUE,
                         admm_local_adapt_niter = 10, ## This spans rho=0.1 to
                                                      ## 100, which is
                                                      ## reasonable.
                         admm_niter = (if(admm_local_adapt)1E3 else 1E4)
                         ## always_first_iter## temporary
                         ){## Basic checks

  ## Basic checks
  if(!is.null(maxdev)){ assertthat::assert_that(maxdev!=0) }
  ## assert_that(!(is.data.frame(ylist[[1]])))
  assertthat::assert_that(sum(is.na(X)) == 0)
  assertthat::assert_that(length(ylist) == nrow(X))
  ## assertthat::assert_that(pie_lambda > 0)
  if(ridge) assert_that(!admm) ## temporary
  assertthat::assert_that(numclust > 1)

  ## Setup
  TT = length(ylist)
  dimdat = ncol(ylist[[1]])
  p = ncol(X)
  if(is.null(mn)) mn = init_mn(ylist, numclust, TT, dimdat, countslist)
  ntlist = sapply(ylist, nrow)
  N = sum(ntlist)

  ## Initialize some objects
  pie = calc_pie(TT, numclust) ## Let's just say it is all 1/K for now.
  denslist_by_clust <- NULL
  objectives = c(+1E20, rep(NA, niter-1))
  sigma = init_sigma(ylist, numclust, TT, fac=sigma.fac) ## (T x numclust x dimdat x dimdat)
  sigma_eig_by_clust = NULL
  zero.betas = zero.alphas = list()
  admm_niters = list()

  ## Warm startable variables
  betas = NULL
  Zs = NULL
  wvecs = NULL
  uws = NULL
  Uzs = NULL

  ## The least elegant solution I can think of.. used only for blocked cv
  if(!is.null(countslist_overwrite)) countslist = countslist_overwrite
  if(!is.null(countslist)) check_trim(ylist, countslist)

  start.time = Sys.time()
  for(iter in 2:niter){
    if(verbose){
      printprogress(iter-1, niter-1, "EM iterations.", start.time = start.time)
    }
    resp <- Estep_covar(mn, sigma, pie, ylist = ylist, numclust = numclust,
                        denslist_by_clust = denslist_by_clust,
                        first_iter = (iter == 2), countslist = countslist)

    ## M step (three parts)
    ## 1. Alpha
    res.alpha = Mstep_alpha(resp, X, numclust, lambda = pie_lambda,
                            zerothresh = zerothresh)
                            ## iter=iter)
    pie = res.alpha$pie
    alpha = res.alpha$alpha
    rm(res.alpha)

    ## print(iter)
    ## if(iter==2){
    ##   save(resp, ylist, X, sigma, numclust, maxdev, sigma_eig_by_clust,
    ##        ## file = file.path("~/Desktop/test-admm-large-T-numclust-10.Rdata"))
    ##        file = file.path("~/Desktop/test-admm-highdim.Rdata"))
    ##   return()
    ## }
    ## load(file.path("~/Desktop/test-admm.Rdata"))

    ## 2. Beta
    if(admm){
      first_iter = (iter == 2)

      ## if(!always_first_iter) first_iter = (iter == 2)
      ## if(always_first_iter) first_iter=TRUE
      res.beta = Mstep_beta_admm(resp, ylist, X,
                                 mean_lambda = mean_lambda,
                                 ## first_iter = (iter == 2),
                                 ## first_iter=TRUE,
                                 first_iter = first_iter,
                                 ## em_iter = iter,
                                 sigma_eig_by_clust = sigma_eig_by_clust,
                                 sigma = sigma, maxdev = maxdev, rho = admm_rho,

                                 betas = betas,
                                 Zs = Zs,
                                 wvecs=wvecs,
                                 uws=uws,
                                 Uzs=Uzs,

                                 err_rel = admm_err_rel,
                                 err_abs = admm_err_abs,
                                 niter = admm_niter,
                                 local_adapt = admm_local_adapt,
                                 local_adapt_niter = admm_local_adapt_niter)
      admm_niters[[iter]] = unlist(res.beta$admm_niters)
    } else {
      res.beta = Mstep_beta(resp, ylist, X,
                             mean_lambda = mean_lambda,
                             sigma = sigma,
                             maxdev = maxdev,
                             ## Temporary
                             ridge = ridge,
                             ridge_lambda = ridge_lambda,
                            ## End of temporary
                             sigma_eig_by_clust = sigma_eig_by_clust,
                             first_iter = (iter == 2),
                             cvxr_ecos_thresh = mstep_cvxr_ecos_thresh,
                             cvxr_scs_eps = mstep_cvxr_scs_eps,
                             zerothresh = zerothresh)
    }

    ## Harvest means
    mn = res.beta$mns
    betas = beta = res.beta$betas

    ## Harvest other things for next iteration's ADMM.
    Zs = res.beta$Zs
    wvecs = res.beta$wvecs
    uws = res.beta$uws
    Uzs = res.beta$Uzs
    rm(res.beta)

    ## Check if the number of zeros in the alphas and betas have stabilized.
    zero.betas[[iter]] = lapply(beta, function(mybeta) which(mybeta==0))
    zero.alphas[[iter]] = which(alpha==0)
    sym_diff <- function(a,b) unique(c(setdiff(a,b), setdiff(b,a)))
    if(zero_stabilize & iter >= 30){ ## If 5 is to low, try 10 instead of 5.

      ## Temporary print message, to see sparsity.
      cat(fill = TRUE)
      print('#Zeros in beta is')
      for( b in beta){
        cat(sum(b[-1,]!=0), "out of", length(b[-1,]), fill=TRUE)
      }
      print('And #zero in alpha is')
      cat(sum(alpha[,-1]!=0), "out of", length(alpha[,-1]), fill=TRUE)
      ## End of temporary

      beta.sym.diffs = Map(sym_diff, zero.betas[[iter]], zero.betas[[iter-1]])
      ## sym_diff(zero.betas[[iter]][[3]], zero.betas[[iter-1]][[3]])
      num.beta.sym.diffs = sapply(beta.sym.diffs, length)
      zero.beta.stable = all(num.beta.sym.diffs <= 1)
      zero.alpha.stable = (length(sym_diff(zero.alphas[[iter]], zero.alphas[[iter-1]])) <= 1)
      if(zero.alpha.stable & zero.beta.stable) break
    }

    ## 3. Sigma
    sigma = Mstep_sigma_covar(resp,
                              ylist,
                              mn,
                              numclust)

    ## 3. (Continue) Decompose the sigmas.
    sigma_eig_by_clust <- eigendecomp_sigma_array(sigma)
    denslist_by_clust <- make_denslist_eigen(ylist, mn, TT, dimdat, numclust,
                                             sigma_eig_by_clust,
                                             countslist)

    ## Calculate the objectives
    objectives[iter] = objective(mn, pie, sigma, ylist,
                                 pie_lambda = pie_lambda,
                                 mean_lambda = mean_lambda,
                                 alpha = alpha, beta = beta,
                                 denslist_by_clust = denslist_by_clust,
                                 countslist = countslist)

    ## print(gc())

    ########################
    ## Make plots ##########
    ########################
    if(plot){
      plot_iter(ylist, countslist, iter, tt=1, TT, mn, sigma, pie, objectives,
                saveplot = TRUE,
                ## saveplot = FALSE,
                plotdir = plotdir)
    }

    ## Check convergence
    ## if(iter > 10){ ## don't stop super early. ## We might not need this.
      if(check_converge_rel(objectives[iter-1],
                            objectives[iter],
                            tol = tol_em)) break
    ## }
    ## if(objectives[iter] > objectives[iter-1] * 1.01 ) break # Additional stopping
                                        ## of the likelihood
                                        ## increasing more
                                        ## than 1%.
  }

  ## Measure time
  lapsetime = difftime(Sys.time(), start.time, units = "secs")
  time_per_iter = lapsetime / (iter-1)

  ## Also calculate per-cytogram likelihoods (NOT divided by nt)
  loglikelihoods = objective(mn, pie, sigma, ylist,
                             pie_lambda = pie_lambda,
                             mean_lambda = mean_lambda,
                             alpha = alpha, beta = beta,
                             denslist_by_clust = denslist_by_clust,
                             countslist = countslist,
                             each = TRUE)

  ## Also reformat the coefficients
  obj <- reformat_coef(alpha, beta, p, numclust, dimdat, X)
  alpha = obj$alpha
  beta = obj$beta

  return(structure(list(alpha = alpha,
                        beta = beta,
                        mn = mn,
                        pie = pie,
                        sigma = sigma,
                        ## denslist_by_clust = denslist_by_clust,
                        objectives = objectives[2:iter],
                        final.iter = iter,
                        time_per_iter = time_per_iter,
                        total_time = lapsetime,
                        loglikelihoods = loglikelihoods,
                        ## Above is output, below are data/algorithm settings.
                        dimdat = dimdat,
                        TT = TT,
                        N = N,
                        p = p,
                        numclust = numclust,
                        X = X,
                        pie_lambda = pie_lambda,
                        mean_lambda = mean_lambda,
                        maxdev=maxdev,
                        refit = refit,
                        niter = niter,
                        admm_niters = admm_niters
                        ), class = "covarem"))
}

## ## Some tests to add
## ## object is the result of having run covarem() or covarem_once().
## check_size <- function(obj){
##   assert_that(check_beta_size(res$beta, p, dimdat, numclust))
##   assert_that(check_alpha_size(res$alpha, p, dimdat))
## }

## check_beta_size <- function(beta, p, dimdat, numclust){
##   all.equal(dim(beta), c(p+1, dimdat, numclust))
## }
## check_alpha_size <- function(alpha, p, dimdat){
##   all.equal(dim(alpha), c(dimdat, p+1))
## }



##' Prediction: given  new X's,  generate a set of means and pies (and return
##' the same Sigma)
##' @param res object returned from covariate EM covarem().
##' @param newx New covariate.
##' @return List containing mean, pie, and sigma.
predict.covarem <- function(res, newx = NULL){

  ## ## Check the dimensions
  ## stopifnot(ncol(new.x) == ncol(res$X))
  ## newx = X[1,,drop=FALSE]
  if(is.null(newx)){
    newx = res$X
  }

  ## Augment it with a dummy variable 1
  if(nrow(newx)>1){
    newx.a = cbind(rep(1, nrow(newx)), newx)
  } else {
    newx.a = c(1, newx)
  }

  TT = nrow(newx) ## This used to be nrow(X)..
  numclust = res$numclust
  dimdat = res$dimdat

  ## Predict the means (manually).
  newmn = lapply(1:numclust, function(iclust){
    newx.a  %*%  res$beta[[iclust]]
  })
  newmn_array = array(NA, dim=c(TT, dimdat, numclust))
  for(iclust in 1:numclust){ newmn_array[,,iclust] = newmn[[iclust]] }

  ## Predict the pies.
  ## newpie = predict(res$alpha.fit, newx=newx, type='response')[,,1]

  piehatmat = as.matrix(exp(cbind(1,newx) %*% t(res$alpha)))
  newpie = piehatmat / rowSums(piehatmat)
  ## predict(fit, newx=X, type="response")[,,1]
  stopifnot(all(dim(newpie) == c(TT,numclust)))
  stopifnot(all(newpie >= 0))

  ## Return all three things
  return(list(newmn = newmn_array,
              newpie = newpie,
              sigma = res$sigma,
              TT = res$TT,
              N = res$N))
}



##' Helper for making list of densities. Returns list by cluster then time
##' e.g. access by \code{denslist_by_clust[[iclust]][[tt]]}
##' @param ylist T-length list each containing response matrices of size (nt x
##'   3), which contains coordinates of the 3-variate particles, organized over
##'   time (T) and with (nt) particles at every time.
##' @param mu (T x dimdat x numclust) array.
##' @param dimdat dimension of data.
##' @param numclust number of clusters.
##' @param TT number of time points
##' @param sigma_eig_by_clust Result of running
##'   \code{eigendecomp_sigma_array(sigma.list[[iter]])}.
##' @return numclust-lengthed list of TT-lengthed.
make_denslist_eigen <- function(ylist, mu,
                                TT, dimdat, numclust,
                                sigma_eig_by_clust,
                                countslist){ ## Temporary

  ## Basic checks
  assert_that(!is.null(sigma_eig_by_clust))

  ## Calculate densities (note to self: nested for loop poses no problems)
  lapply(1:numclust, function(iclust){
    mysigma_eig <- sigma_eig_by_clust[[iclust]]
      lapply(1:TT, function(tt){
        return(dmvnorm_fast(ylist[[tt]],
                            mu[tt,,iclust],
                            sigma_eig=mysigma_eig))
    })
  })
}

## ## helper function for plotting within covarem_once() iterations.
## plot_iter_3d <- function(ylist, countslist, iter=NULL, tt = 1, TT, mn, sigma, pie, saveplot=TRUE){

##   if(!is.null(countslist)){
##     cex = (countslist[[1]]/max(countslist[[1]]))*5+.5
##   } else {
##     cex = .5
##   }
##   ylist_collapsed = do.call(rbind, ylist)
##   ntsum = nrow(ylist_collapsed)
##   ylist_collapsed = ylist_collapsed[sample(1:ntsum, 10000),]
##   if(saveplot) png(file=file.path(plotdir,
##                      paste0("iteration-", iter, ".png")), width=1200, height=1200)
##   par(mfrow=c(2,2))
##   dimslist = list(1:2,2:3,c(3,1))
##   for(dims in dimslist){
##     names = c("fsc_small", "chl_small","pe")
##     plot(x=ylist[[tt]][,dims[1]],
##          y=ylist[[tt]][,dims[2]],
##          ## col='grey50',
##          col=rgb(0,0,0,0.1),
##          pch=16,
##          cex = cex,
##          ## col='pink',
##          type='p',
##          xlab = names[dims[1]],
##          ylab = names[dims[2]],
##          cex.lab = 1.5,
##          cex.axis = 1.5
##          )

##     ## this time point's mean
##     cols = 1:numclust
##     points(x = mn[tt,dims[1],],
##            y = mn[tt,dims[2],],
##            pch = 4, col = cols,
##            ## cex = pie[1,]/max(pie[1,])*10)
##            cex = log(pie[1,]/max(pie[1,]) + 1)*3)

##     ## All time points' means
##     for(tt in 1:TT){
##       for(kk in 1:numclust){
##         points(x=mn[tt,dims[1],kk],
##                y=mn[tt,dims[2],kk],
##                pch=16, col=cols[kk],
##                cex = 1)
##       }
##     }

##     for(kk in 1:numclust){
##       lines(ellipse::ellipse(x = sigma[kk,dims,dims],
##                              centre = mn[tt,dims, kk]),
##             lwd=1, col=cols[kk], lty=1)
##     }
##   }
##   plot(objectives[2:(min(100, length(objectives)))],
##        type= 'o',
##        ylab = "Objective Value", xlab = "EM Iteration",
##        cex.lab = 1.5,
##        cex.axis = 1.5)
##   if(saveplot)  graphics.off()
## }


## helper function for plotting within covarem_once() iterations.
plot_iter <- function(ylist, countslist, iter=NULL, tt = 1, TT, mn, sigma, pie, objectives,
                      saveplot=TRUE, plotdir=NULL){

  if(!is.null(countslist)){
    cex = (countslist[[tt]]/max(countslist[[tt]]))*5+.5
  } else {
    cex = .5
  }

  if(saveplot){
    png(file=file.path(plotdir,
                     paste0("iteration-", iter, ".png")), width=1200, height=1200)
    print(  file.path(plotdir,
          paste0("iteration-", iter, ".png")))
  }

  dimdat = ncol(ylist[[1]])
  if(dimdat==3){
    par(mfrow=c(2,2))
    dimslist = list(1:2, 2:3, c(3,1))
  }
  if(dimdat==2){
    par(mfrow=c(1,2))
    dimslist = list(1:2)
  }

  ylist_collapsed = do.call(rbind, ylist)
  xlims = lapply(dimslist, function(dims) range(ylist_collapsed[,dims[1]]))
  ylims = lapply(dimslist, function(dims) range(ylist_collapsed[,dims[2]]))
  ntsum = nrow(ylist_collapsed)
  if(ntsum>10000) ylist_collapsed = ylist_collapsed[sample(1:ntsum, 10000),]

  ## Make the plots
  for(ii in 1:length(dimslist)){
    dims = dimslist[[ii]]
    xlim = xlims[[ii]]
    ylim = ylims[[ii]]
    names = c("fsc_small", "chl_small","pe")
    plot(x=ylist[[tt]][,dims[1]],
         y=ylist[[tt]][,dims[2]],
         ## col='grey50',
         col=rgb(0,0,0,0.1),
         pch=16,
         cex = cex,
         ## col='pink',
         type='p',
         xlab = names[dims[1]],
         ylab = names[dims[2]],
         cex.lab = 1.5,
         cex.axis = 1.5,
         xlim = xlim,
         ylim = ylim
         )

    ## this time point's mean
    cols = 1:numclust
    points(x = mn[tt,dims[1],],
           y = mn[tt,dims[2],],
           pch = 16, col = cols,
           cex = pie[1,]/max(pie[1,])*5)
           ## cex = log(pie[1,]/max(pie[1,]) + 1)*3)


    ## All time points' means
    for(ttt in 1:TT){
      for(kk in 1:numclust){
        points(x=mn[ttt,dims[1],kk],
               y=mn[ttt,dims[2],kk],
               pch=16, col=cols[kk],
               cex = .2)
      }
    }

    ## Draw the ellipses
    for(kk in 1:numclust){
      lines(ellipse::ellipse(x = sigma[kk,dims,dims],
                             centre = mn[tt,dims, kk]),
            lwd=1, col=cols[kk], lty=1)
    }
  }
  plot(objectives[2:(min(100, length(objectives)))],
       type = 'o',
       ylab = "Objective Value", xlab = "EM Iteration",
       cex.lab = 1.5,
       cex.axis = 1.5)
  if(saveplot)  graphics.off()
}
