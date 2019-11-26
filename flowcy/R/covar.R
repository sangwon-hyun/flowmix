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
##' @return List containing fitted parameters and means and mixture weights,
##'   across algorithm iterations. beta is a (p+1 x 3 x dimdat) array. Alpha is
##'   a (dimdat x (p+1)) array.
##'
##' @export
covarem_once <- function(ylist, X,
                         ## ylist_orig,## temporary
                         countslist = NULL,
                         numclust, niter = 100,
                         mn = NULL, pie_lambda = 0,
                         mean_lambda = 0, verbose = FALSE,
                         warmstart  =  c("none", "rough"), sigma.fac = 1, tol = 1E-3,
                         refit = FALSE, ## EXPERIMENTAL FEATURE.
                         sel_coef = NULL,
                         maxdev = NULL,
                         manual.bin = FALSE,
                         manual.grid = NULL,
                         countslist_overwrite = NULL
                         ## ridge = FALSE,
                         ## ridge_lambda = 0
                         ){
  ## Basic checks
  if(!is.null(maxdev)){ assert_that(maxdev!=0) } ## Preventing the maxdev=FALSE mistake.
  ## assert_that(!(is.data.frame(ylist[[1]])))

  ## Setup.
  TT = length(ylist)
  dimdat = ncol(ylist[[1]])
  p = ncol(X)
  if(is.null(mn)) mn = init_mn(ylist, numclust, TT, dimdat, warmstart)
  if(!is.null(mn)) orig.mn = mn
  ntlist = sapply(ylist, nrow)
  bin = !is.null(countslist)

  ## Initialize some objects
  pie = calc_pie(TT, numclust) ## Let's just say it is all 1/K for now.
  denslist_by_clust <- NULL
  objectives = c(+1E20, rep(NA, niter-1))
  sigma = init_sigma(ylist, numclust, TT, fac=sigma.fac) ## (T x numclust x dimdat x dimdat)
  sigma_eig_by_clust = NULL

  ## If binning manual,
  if(manual.bin){
    assert_that(!is.null(manual.grid))
    bin = TRUE

    ## Bin data
    cat("binning started", fill=TRUE)
    print(Sys.time())
    ## reslist = lapply(ylist, bin_one_cytogram, manual.grid)
    start.time1 = Sys.time()
    reslist = lapply(1:TT, function(tt){
      if(tt %% 100==0)printprogress(tt, TT, "binning", start.time = start.time1, fill=TRUE)
      bin_one_cytogram(ylist[[tt]], manual.grid)})
    ylist = lapply(reslist, function(res) res$ybin)
    countslist = lapply(reslist, function(res) res$counts)
    print(Sys.time())
    cat("binning ended", fill=TRUE)
  }
  ## Endof alternative

  if(!is.null(countslist_overwrite))countslist = countslist_overwrite ## The least elegant solution I can think of..


  ## Temporary: Also, if we are using binned counts, make sure that ylist and
  ## countslist are trimmed i.e. make sure that in
  if(bin) check_trim(ylist, countslist)

  start.time = Sys.time()
  for(iter in 2:niter){
    if(verbose) printprogress(iter-1, niter-1, "EM iterations.", start.time = start.time)

    start_time_per_iter = Sys.time()

    resp <- Estep_covar(mn, sigma, pie,
                        ylist = ylist,
                        numclust,
                        denslist_by_clust = denslist_by_clust,
                        first_iter = (iter == 2),
                        countslist = countslist)

    ## If countslist is provided, further weight them.
    if(!is.null(countslist)){
      resp <- Map(function(myresp, mycount){
        myresp * mycount
      }, resp, countslist)
    }

    ## Conduct M step
    ## 1. Alpha
    res.alpha = Mstep_alpha(resp,
                            X, numclust,
                            lambda = pie_lambda,
                            refit = refit,
                            sel_coef = sel_coef,
                            bin = bin)
    pie = res.alpha$pie
    alpha = res.alpha$alpha
    rm(res.alpha)

    ## 2. Beta
    res.beta = Mstep_beta(resp, ylist, X,
                          mean_lambda = mean_lambda,
                          sigma,
                          refit = refit,
                          sel_coef = sel_coef, maxdev = maxdev,
                          sigma_eig_by_clust = sigma_eig_by_clust,
                          first_iter = (iter == 2),
                          ## ridge = ridge,
                          ## ridge_lambda = ridge_lambda,
                          ## ridge_pie = pie,
                          bin = bin)
    mn = res.beta$mns
    beta = res.beta$beta
    rm(res.beta)

    ## 3. Sigma
    sigma = Mstep_sigma_covar(resp,
                              ylist,
                              mn,
                              numclust,
                              bin = bin)

    ## 3. (Continue) Decompose the sigmas.
    sigma_eig_by_clust <- eigendecomp_sigma_array(sigma)
    denslist_by_clust <- make_denslist_eigen(ylist,
                                             mn,
                                             TT, dimdat,
                                             numclust, sigma_eig_by_clust,
                                             bin=bin,
                                             countslist)

    ## Calculate the objectives
    objectives[iter] = objective_overall_cov(mn,
                                             pie,
                                             sigma,
                                             TT,
                                             dimdat,
                                             numclust,
                                             ylist,
                                             pie_lambda = pie_lambda,
                                             mean_lambda = mean_lambda,
                                             alpha = alpha,
                                             beta = beta,
                                             denslist_by_clust = denslist_by_clust,
                                             countslist = countslist)

    ########################
    ## Make plots ##########
    ########################
    ## browser()
    par(mfrow=c(1,2))
    if(!is.null(countslist)){
      cex = (countslist[[1]]/max(countslist[[1]]))*5+.5
    } else {
      cex = .5
    }
    ylist_collapsed = do.call(rbind, ylist)
    ntsum = nrow(ylist_collapsed)
    ylist_collapsed = ylist_collapsed[sample(1:ntsum, 10000),]
    plot(x = ylist_collapsed[,1],
         y = ylist_collapsed[,2],
         pch=16,
         ## col='grey80',
         col=rgb(0,0,0,0.1),
         cex = cex,
         type='p')
    tt = 1
    points(x=ylist[[tt]][,1],
           y=ylist[[tt]][,2], pch=16,
           col='grey50',
            cex = cex,
            type='p')

    ## for(iclust in 1:numclust){
    ##   x = mn[,1, iclust]
    ##   y = mn[,2, iclust]
    ##   points(x=as.numeric(x), y=as.numeric(y),col=iclust)
    ## }
    cols = 1:numclust
    points(x=mn[tt,1,],
           y=mn[tt,2,],
           pch="x", col=cols,
           cex = pie[1,]/max(pie[1,])*10)

    for(tt in 1:TT){
    points(x=mn[tt,1,],
           y=mn[tt,2,],
           pch=16, col="blue",
           cex = .3)
    points(x=orig.mn[tt,1,],
           y=orig.mn[tt,2,],
           pch=16, col="red",
           cex = .3)
    }

    for(kk in 1:numclust){
      lines(ellipse::ellipse(x = sigma[kk,,],
                             centre = mn[tt,, kk]
                             ), lwd=1, col=cols[kk], lty=1)
    }
    plot(objectives[2:50], type='o') ##temporary
    ## End of plotting  ##########

    ## Check convergence
    if(check_converge_rel(objectives[iter-1],
                          objectives[iter],
                          tol = tol)) break
  }

  ## Threshold (now that we're using CVXR for beta)
  betathresh = 1E-3
  beta = lapply(beta, function(a){
    a[abs(a)<betathresh] = 0
    a
  })

  return(structure(list(alpha = alpha,
                        beta = beta,
                        mn = mn,
                        pie = pie,
                        sigma = sigma,
                        resp = resp, ## Temporary
                        denslist_by_clust = denslist_by_clust, ## Temporary
                        objectives = objectives[1:iter],
                        final.iter = iter,
                        ## Above is output, below are data/algorithm settings.
                        dimdat = dimdat,
                        TT = TT,
                        p = p,
                        numclust = numclust,
                        X = X,
                        pie_lambda = pie_lambda,
                        mean_lambda = mean_lambda,
                        maxdev=maxdev,
                        refit=refit,
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
