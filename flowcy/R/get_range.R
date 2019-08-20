##' Getting the maximum lambda values from the first EM iteration. Taken
##' directly from main function for covariate EM, as of June 23rd 2019.
##' @param ylist T-length list each containing response matrices of size (nt x
##'   3), which contains coordinates of the 3-variate particles, organized over
##'   time (T) and with (nt) particles at every time.
##' @param X Matrix of size (T x p+1)
##' @param pie.list (T by K)
##' @param mean_lambda lambda for lasso for the mean.
##' @param pie_lambda lambda for lasso for pie.
##' @return List containing fitted parameters and means and mixture weights,
##'   across algorithm iterations.
covarem_getrange <- function(ylist, X=NULL, numclust, niter=100, mn=NULL, pie_lambda=0,
                             mean_lambda=0, verbose=FALSE,
                             warmstart = c("none", "rough"), sigma.fac=1, tol=1E-6){

  ## Setup.
  ntlist = sapply(ylist, nrow)
  dimdat = ncol(ylist[[1]])
  TT = length(ylist)
  p = ncol(X)
  warmstart = match.arg(warmstart)

  ## Initialize.
  beta = init_beta(TT, p, dimdat, numclust)
  alpha = init_alpha(dimdat, p)
  if(is.null(mn)){
    if(warmstart=="rough"){
      mn = warmstart_covar(ylist, numclust)
    } else if (warmstart=="none"){
      mn = init_mu(lapply(ylist, cbind), numclust, TT)
    } else {
      stop("warmstart option not recognized")
    }
  }
  pie = calc_pie(TT, numclust) ## Let's just say it is all 1/K for now.
  sigma = init_sigma(ylist, numclust, TT, fac=sigma.fac) ## (T x numclust x dimdat x dimdat)

  ## Initialize alpha and beta
  beta.list = alpha.list = sigma.list = pie.list = mn.list = list()
  objectives = rep(NA, niter)
  objectives[1] = -1E20 ## Fake
  beta.list[[1]] = beta ## beta.list: Each element is a (T x p+1 x 3 x K) array
  alpha.list[[1]] = alpha ## alpha.list: Each element is a T by p+1 array
  mn.list[[1]] = mn
  sigma.list[[1]] = sigma
  pie.list[[1]] = pie

  start.time=Sys.time()
  for(iter in 2){
    if(verbose) printprogress(iter, niter, "EM iterations.", start.time=start.time)

    ## Conduct E step
    resp <- Estep_covar(mn.list[[iter-1]],
                        sigma.list[[iter-1]],
                        pie.list[[iter-1]],
                        ylist,
                        numclust)  ## This should be (T x numclust x dimdat x dimdat)

    ## Conduct M step
    ## 1. Alpha
    max_lambda_alpha = Mstep_alpha_getrange(resp,
                                            X, numclust,
                                            lambda=pie_lambda)

    ## 2. Beta
    max_lambda_beta = Mstep_beta_getrange(resp, ylist, X,
                                          mean_lambda=mean_lambda,
                                          sigma.list[[iter-1]])

  }
  return(list(max_lambda_beta=max_lambda_beta,
              max_lambda_alpha=max_lambda_alpha))

}



##' Getting the range of lambda values from the first alpha update in the
##' covariate EM algorithm, using GLMnet as the workhorse.
##' @param resp Is an (T x nt x K) array.
##' @param X Covariate matrix (T x dimdat).
##' @return The multinomial logit model coefficients. A matrix of dimension (K x
##'   (p+1)).
Mstep_alpha_getrange <- function(resp, X, numclust, lambda=0, alpha=1){

  TT = nrow(X)
  p = ncol(X)

  ## Calculate the summed responsibilities
  resp.sum = t(sapply(resp, colSums)) ## (T x numclust)

  stopifnot(dim(resp)==c(TT, numclust))

  ## Fit the model, possibly with some regularization
  fit = glmnet::glmnet(x=X, y=resp.sum, family="multinomial",
                       alpha=alpha)
  return(max(fit$lambda))
}



##' Getting the range of lambda values from the first M step of beta using a faster
##' lasso regression formulation.
Mstep_beta_getrange <- function(resp, ylist, X, mean_lambda=0, sigma, numclust){

  ## Preliminaries
  TT = length(ylist)
  numclust = ncol(resp[[1]])
  dimdat = ncol(ylist[[1]])
  ntlist = sapply(ylist, nrow)

  Xa = cbind(rep(1, TT), X)

  ## Setup
  manip_obj = manip(ylist, X, resp, sigma, numclust)
  Xtildes = manip_obj$Xtildes
  yvecs = manip_obj$yvecs

  ## Intercepts are to be excluded.
  exclude.from.penalty = (0:(dimdat-1))*(ncol(Xa)) + 1

  ## Give the glmnet function pre-calculated Y's and X's.
  max_lambdas_by_clust = lapply(1:numclust, function(iclust){

    ## Give the glmnet function pre-calculated Y's and X's.
    fit = glmnet::glmnet(x=Xtildes[[iclust]], y = yvecs[[iclust]],
                         alpha=1, intercept = FALSE, family = "gaussian")
    max_lambda = max(fit$lambda)
    return(max_lambda)
  })
  return(max_lambdas_by_clust)
}



##' Estimate maximum lambda values.
##' @param ylist List of responses.
##' @param X Covariates.
##' @param numclust Number of clusters.
##' @param maxfac Defaults to 32.
##' @param ... Other arguments to \code{covarem_once()}.
##' @return list containing the two maximum values to use.
get_max_lambda <- function(ylist, X, numclust, maxfac=32, ...){

  ## Get range of regularization parameters.
  res0 = covarem_getrange(ylist=ylist, X=X, numclust=numclust, niter=2)

  fac = 2
  while(fac <= maxfac){

    ## Checking the maximum lambda value
    max_lambda_beta = max(unlist(res0$max_lambda_beta)) * fac
    max_lambda_alpha = max(res0$max_lambda_alpha) * fac

    ## Checking that the max actually zeros out.
    res = covarem(ylist=ylist, X=X, numclust=numclust,
                  mean_lambda=max_lambda_beta, pie_lambda=max_lambda_alpha, ...)
    alpha.checks.out = all(res$alpha[,-1]==0)
    beta.checks.out = all(sapply(res$beta, function(cf){ all(cf[-1,]==0) }))
    if(alpha.checks.out & beta.checks.out) break
    fac = fac * 2
  }
  return(list(beta=max_lambda_beta, alpha=max_lambda_alpha))
}
