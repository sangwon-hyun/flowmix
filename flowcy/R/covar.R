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
##' @param tol_em Relative tolerance for EM convergence.
##' @return List containing fitted parameters and means and mixture weights,
##'   across algorithm iterations. beta is a (p+1 x 3 x dimdat) array. Alpha is
##'   a (dimdat x (p+1)) array.
##'
##' @export
covarem_once <- function(ylist, X,
                         countslist = NULL,
                         numclust, niter = 1000,
                         mn = NULL, pie_lambda = 0,
                         mean_lambda = 0, verbose = FALSE,
                         sigma.fac = 1, tol_em = 1E-5,
                         refit = FALSE, ## EXPERIMENTAL FEATURE.
                         sel_coef = NULL,
                         maxdev = NULL,
                         manual.bin = FALSE,
                         manual.grid = NULL,
                         countslist_overwrite = NULL,
                         zero_stabilize  = FALSE,
                         ## ridge = FALSE,
                         ## ridge_lambda = 0
                         plot = FALSE,
                         plotdir = "~/Desktop",
                         init_mn_flatten = FALSE,
                         ## beta Mstep (CVXR) settings
                         mstep_cvxr_ecos_thresh = 1E-8,
                         mstep_cvxr_scs_eps = 1E-5,
                         zerothresh = 1E-6,
                         ## beta Mstep (ADMM) settings
                         admm = TRUE,
                         admm_rho = 10,
                         admm_err_rel = 1E-3,
                         ## beta M step (Locally Adaptive ADMM) settings
                         admm_local_adapt = TRUE,
                         admm_local_adapt_niter = 5,
                         admm_niter = (if(admm_local_adapt)3E2 else 1E4)
                         ){## Basic checks

  if(!is.null(maxdev)){ assert_that(maxdev!=0) } ## Preventing the maxdev=FALSE mistake.
  ## assert_that(!(is.data.frame(ylist[[1]])))
  assert_that(sum(is.na(X)) == 0)
  assert_that(length(ylist) == nrow(X))


  ## Setup.
  TT = length(ylist)
  dimdat = ncol(ylist[[1]])
  p = ncol(X)
  if(is.null(mn)) mn = init_mn(ylist, numclust, TT, dimdat, countslist)
  ntlist = sapply(ylist, nrow)
  bin = !is.null(countslist)

  ## Initialize some objects
  pie = calc_pie(TT, numclust) ## Let's just say it is all 1/K for now.
  denslist_by_clust <- NULL
  objectives = c(+1E20, rep(NA, niter-1))
  sigma = init_sigma(ylist, numclust, TT, fac=sigma.fac) ## (T x numclust x dimdat x dimdat)
  sigma_eig_by_clust = NULL
  zero.betas = zero.alphas = list()

  ## The least elegant solution i can think of.. used for blocked cv
  if(!is.null(countslist_overwrite))countslist = countslist_overwrite
  if(bin) check_trim(ylist, countslist)

  start.time = Sys.time()
  for(iter in 2:niter){
    if(verbose){
      printprogress(iter-1, niter-1, "EM iterations.", start.time = start.time)
    }

    resp <- Estep_covar(mn, sigma, pie, ylist = ylist, numclust = numclust,
                        denslist_by_clust = denslist_by_clust,
                        first_iter = (iter == 2), countslist = countslist)

    ## If countslist is provided, further weight them.
    if(!is.null(countslist)){
      resp <- Map(function(myresp, mycount){ myresp * mycount },
                  resp, countslist)
    }

    ## M step (three parts)
    ## 1. Alpha
    res.alpha = Mstep_alpha(resp, X, numclust, lambda = pie_lambda,
                            refit = refit, sel_coef = sel_coef, bin = bin,
                            ## thresh = thresh,
                            cvxr_ecos_thresh = mstep_cvxr_ecos_thresh,
                            zerothresh = zerothresh)
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
      res.beta = Mstep_beta_admm(resp, ylist, X,
                                 mean_lambda = mean_lambda,
                                 first_iter = (iter == 2),
                                 sigma_eig_by_clust = sigma_eig_by_clust,
                                 sigma = sigma, maxdev = maxdev, rho = admm_rho,
                                 err_rel = admm_err_rel,
                                 niter = admm_niter,
                                 local_adapt = admm_local_adapt)
    } else {
      res.beta = Mstep_beta(resp, ylist, X,
                            mean_lambda = mean_lambda,
                            sigma = sigma,
                            maxdev = maxdev,
                            ## This is fluff:
                            ## refit = refit,
                            ## sel_coef = sel_coef,
                            ## end of fluff
                            sigma_eig_by_clust = sigma_eig_by_clust,
                            first_iter = (iter == 2), bin = bin,
                            cvxr_ecos_thresh = mstep_cvxr_ecos_thresh,
                            cvxr_scs_eps = mstep_cvxr_scs_eps,
                            zerothresh = zerothresh)
    }

    mn = res.beta$mns
    beta = res.beta$beta
    rm(res.beta)

    ## Check if the number of zeros in the alphas and betas have stabilized.
    zero.betas[[iter]] = lapply(beta, function(mybeta) which(mybeta==0))
    zero.alphas[[iter]] = which(alpha==0)
    sym_diff <- function(a,b) unique(c(setdiff(a,b), setdiff(b,a)))
    if(zero_stabilize & iter >= 5){ ## If 5 is to low, try 10 instead of 5.

      ## Temporary print message, to see sparsity.
      cat(fill = TRUE)
      print('beta')
      for( b in beta){
        cat(sum(b[-1,]!=0), "out of", length(b[-1,]), fill=TRUE)
      }
      print('alpha')
      cat(sum(alpha[,-1]!=0), "out of", length(alpha[,-1]), fill=TRUE)
      ## End of temporary

      beta.sym.diffs = Map(sym_diff, zero.betas[[iter]], zero.betas[[iter-1]])
      sym_diff(zero.betas[[iter]][[3]], zero.betas[[iter-1]][[3]])
      num.beta.sym.diffs = sapply(beta.sym.diffs, length)
      zero.beta.stable = all(num.beta.sym.diffs <= 1)
      zero.alpha.stable = (length(sym_diff(zero.alphas[[iter]], zero.alphas[[iter-1]])) <= 1)
      if(zero.alpha.stable & zero.beta.stable) break
    }

    ## 3. Sigma
    sigma = Mstep_sigma_covar(resp,
                              ylist,
                              mn,
                              numclust,
                              bin = bin)

    ## 3. (Continue) Decompose the sigmas.
    sigma_eig_by_clust <- eigendecomp_sigma_array(sigma)
    denslist_by_clust <- make_denslist_eigen(ylist, mn, TT, dimdat, numclust,
                                             sigma_eig_by_clust, bin=bin,
                                             countslist)

    ## Calculate the objectives
    objectives[iter] = objective_overall_cov(mn, pie, sigma, TT, dimdat,
                                             numclust, ylist,
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
                plotdir = plotdir)
    }

    ## Check convergence
    if(check_converge_rel(objectives[iter-1],
                          objectives[iter],
                          tol = tol_em)) break
  }

  ## Measure time
  lapsetime = difftime(Sys.time(), start.time, units = "secs")
  time_per_iter = lapsetime / (iter-1)

  return(structure(list(alpha = alpha,
                        beta = beta,
                        mn = mn,
                        pie = pie,
                        sigma = sigma,
                        ## resp = resp, ## Temporary
                        denslist_by_clust = denslist_by_clust, ## Temporary
                        objectives = objectives[1:iter],
                        final.iter = iter,
                        time_per_iter = time_per_iter,
                        total_time = lapsetime,
                        ## Above is output, below are data/algorithm settings.
                        dimdat = dimdat,
                        TT = TT,
                        p = p,
                        numclust = numclust,
                        X = X,
                        pie_lambda = pie_lambda,
                        mean_lambda = mean_lambda,
                        maxdev=maxdev,
                        refit = refit,
                        niter = niter
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
              sigma = res$sigma))
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
                                countslist, ## Temporary
                                bin){ ## Temporary

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
  ylist_collapsed = do.call(rbind, ylist)
  ntsum = nrow(ylist_collapsed)
  ylist_collapsed = ylist_collapsed[sample(1:ntsum, 10000),]
  if(saveplot){
    png(file=file.path(plotdir,
                     paste0("iteration-", iter, ".png")), width=1200, height=1200)
    print(  file.path(plotdir,
          paste0("iteration-", iter, ".png")))
    }
  par(mfrow=c(2,2))
  dimslist = list(1:2, 2:3, c(3,1))
  for(dims in dimslist){
    xlim = range(ylist_collapsed[,dims[1]])
    ylim = range(ylist_collapsed[,dims[2]])
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
