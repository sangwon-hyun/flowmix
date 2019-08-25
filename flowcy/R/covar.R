##' Main function for covariate EM. Does covariate EM with |nrep| restarts (5 by
##' default).
##' @param ... Arguments for \code{covarem_once()}.
##' @param nrep Number of restarts.
##' @return The |covarem| class object that had the best likelihood.
covarem <- function(..., nrep=5){

  ## Don't do many replicates if warmstartable mean exists
  dots <- list(...)
  if(!is.null(dots$mn)) nrep = 1

  ## Do |nrep| replicates.
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
##'   across algorithm iterations.
covarem_once <- function(ylist, X = NULL, numclust, niter = 100, mn = NULL, pie_lambda = 0,
                    mean_lambda = 0, verbose = FALSE,
                    warmstart  =  c("none", "rough"), sigma.fac = 1, tol = 1E-6,
                    refit = FALSE, ## EXPERIMENTAL FEATURE.
                    sel_coef = NULL,
                    maxdev = NULL,
                    faster_mvn=FALSE,
                    eigenspeed=FALSE
                    ){

  ## Setup.
  dimdat = ncol(ylist[[1]])
  TT = length(ylist)
  p = ncol(X)
  warmstart = match.arg(warmstart)
  sigma_eig_by_dim <- NULL
  denslist_by_clust <- NULL

  ## Initialize.
  beta = init_beta(TT, p, dimdat, numclust)
  alpha = init_alpha(dimdat, p)
  if(is.null(mn)){
    if(warmstart == "rough"){
      mn = warmstart_covar(ylist, numclust)
    } else if (warmstart == "none"){
      mn = init_mn(lapply(ylist, cbind), numclust, TT)
    } else {
      stop("warmstart option not recognized")
    }
  }
  pie = calc_pie(TT, numclust) ## Let's just say it is all 1/K for now.
  sigma = init_sigma(ylist, numclust, TT, fac=sigma.fac) ## (T x numclust x dimdat x dimdat)

  ## Initialize alpha and beta
  beta.list = alpha.list = sigma.list = pie.list = mn.list = list()
  objectives = rep(NA, niter)
  objectives[1] = +1E20 ## Fake
  beta.list[[1]] = beta ## beta.list: Each element is a (T x p+1 x 3 x K) array
  alpha.list[[1]] = alpha ## alpha.list: Each element is a T by p+1 array
  mn.list[[1]] = mn  ## (T x dimdat x numclust)
  sigma.list[[1]] = sigma
  pie.list[[1]] = pie

  start.time = Sys.time()
  for(iter in 2:niter){
    if(verbose) printprogress(iter, niter, "EM iterations.", start.time = start.time)

    ## Conduct E step
    resp <- Estep_covar(mn.list[[iter-1]],
                        sigma.list[[iter-1]],
                        pie.list[[iter-1]],
                        ylist,
                        numclust,
                        faster_mvn = faster_mvn,
                        sigma_eig_by_dim = sigma_eig_by_dim,
                        denslist_by_clust=denslist_by_clust
                        )  ## This should be (T x numclust x dimdat x dimdat)

    ## Conduct M step
    ## 1. Alpha
    res.alpha = Mstep_alpha(resp,
                            X, numclust,
                            lambda = pie_lambda,
                            refit = refit,
                            sel_coef = sel_coef)
    alpha.list[[iter]] = res.alpha$alpha
    pie.list[[iter]] = res.alpha$pie

    ## 2. Beta
    res.beta = Mstep_beta(resp, ylist, X, mean_lambda = mean_lambda,
                          sigma.list[[iter-1]], refit = refit,
                          sel_coef = sel_coef, maxdev = maxdev,
                          sigma_eig_by_dim = sigma_eig_by_dim)
    beta.list[[iter]] = res.beta$beta
    mn.list[[iter]]    = res.beta$mns

    ## 3. Sigma
    sigma.list[[iter]] <- Mstep_sigma_covar(resp,
                                            ylist,
                                            mn.list[[iter]],
                                            numclust)

    ## 3. (Continue) Eigendecomp the sigmas.
    if(eigenspeed){
      sigma_eig_by_dim <- eigendecomp_sigma_array(sigma.list[[iter]])
    }

    denslist_by_clust <- make_denslist(ylist, mn.list[[iter]],
                                       sigmalist[[iter]], TT,
                                       dimdat, numclust,
                                       sigma_eig_by_dim)

    ## Calculate the objectives
    objectives[iter] = objective_overall_cov(mn.list[[iter]],
                                             pie.list[[iter]],
                                             sigma.list[[iter]],
                                             TT,
                                             dimdat,
                                             numclust,
                                             ylist,
                                             pie_lambda = pie_lambda,
                                             mean_lambda = mean_lambda,
                                             alpha = res.alpha$alpha,
                                             beta = res.beta$beta,
                                             faster_mvn=faster_mvn,
                                             sigma_eig_by_dim = sigma_eig_by_dim,
                                             denslist_by_clust = denslist_by_clust
                                             )

    ## Check convergence
    if(check_converge_rel(objectives[iter-1],
                          objectives[iter], tol = tol)) break
  }

  ## Threshold (now that we're using CVXR for beta)
  beta = beta.list[[iter]]
  betathresh = 1E-3
  beta = lapply(beta, function(a){
    a[abs(a)<betathresh] = 0
    a
  })

  return(structure(list(alpha = alpha.list[[iter]],
                        ## alpha.fit=res.alpha$fit,
                        beta = beta,
                        mn = mn.list[[iter]],
                        pie = pie.list[[iter]],
                        sigma = sigma.list[[iter]],
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


##' Prediction: given  new X's,  generate a set of means and pies (and return
##' the same Sigma)
##' @param res object returned from covariate EM covarem().
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

  TT = length(ylist)
  numclust = res$numclust

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
##' @param sigma_eig_by_dim Result of running
##'   \code{eigendecomp_sigma_array(sigma.list[[iter]])}.
make_denslist <- function(ylist, mu,
                          sigma, TT, dimdat, numclust,
                          sigma_eig_by_dim){ ## This is experimental.

  assert_that(!is.null(sigma_eig_by_dim))
  lapply(1:numclust, function(iclust){

    mysigma_eig = sigma_eig_by_dim[[iclust]]
    lapply(1:TT,function(tt){

      ## Setup
      mydat = ylist[[tt]]
      mymu = mu[tt,,iclust]

      ## Calculate weighted density
      return(dmvnorm_fast(mydat,
                          mu=mymu,
                          sigma_eig=mysigma_eig))
    })
  })
}
