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
covarem_getrange <- function(ylist, X=NULL, numclust, niter=1000, mn=NULL, pie_lambda=0,
                             mean_lambda=0, verbose=FALSE,
                             warmstart = c("none", "rough"), sigma.fac=1, tol=1E-6){

  ## Setup.
  ntlist = sapply(ylist, nrow)
  dimdat = ncol(ylist[[1]])
  TT = length(ylist)
  p = ncol(X)
  warmstart = match.arg(warmstart)

  ## Initialize.
  beta = init_beta(p, dimdat, numclust)
  alpha = init_alpha(dimdat, p)
  if(is.null(mn)){
    if(warmstart=="rough"){
      mn = warmstart_covar(ylist, numclust)
    } else if (warmstart=="none"){
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
                        numclust,
                        first_iter = (iter == 2)
                        )  ## This should be (T x numclust x dimdat x dimdat)

    ## Conduct M step
    ## 1. Alpha
    max_lambda_alpha = Mstep_alpha_getrange(resp,
                                            X, numclust,
                                            lambda=pie_lambda)

    ## 2. Beta
    max_lambda_beta = Mstep_beta_getrange(resp, ylist, X,
                                          mean_lambda=mean_lambda,
                                          sigma.list[[iter-1]]
                                          )

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
  manip_obj = manip(ylist, Xa, resp, sigma, numclust,
                    first_iter = TRUE)

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
                  mean_lambda=max_lambda_beta,
                  pie_lambda=max_lambda_alpha, ...)
    alpha.checks.out = all(res$alpha[,-1]==0)
    beta.checks.out = all(sapply(res$beta, function(cf){ all(cf[-1,]==0) }))
    if(alpha.checks.out & beta.checks.out) break
    fac = fac * 2
  }

  return(list(beta=max_lambda_beta, alpha=max_lambda_alpha))
}



##' Estimate maximum lambda values.
##' @param ylist List of responses.
##' @param X Covariates.
##' @param numclust Number of clusters.
##' @param max_lambda_beta Defaults to 4000.
##' @param max_lambda_beta Defaults to 1000.
##' @param iimax Maximum value of x for 2^{-x} factors to try.
##' @param ... Other arguments to \code{covarem_once()}.
##' @return list containing the two maximum values to use.
##' @examples
##'
##' obj = generate_data_generic(p=5, TT=500, fac=1, nt=7000, dimdat=3)
##' ylist = obj$ylist
##' X = obj$X
##'
##' ## No parallel:
##' maxres = get_max_lambda_new(ylist, X, numclust, verbose=TRUE,
##'                             nrep = 4,
##'                             ## Function settings
##'                             parallelize = FALSE,
##'                             iimax = 20,
##'                             niter = 1000,
##'                             max_lambda_alpha = 10000,
##'                             tol = 1E-3 ## This doesn't need to be so low here.
##'                             )
##' ## Yes parallel:
##' cl = get_cl(3)
##' parallel::clusterExport(cl, ls())
##' parallel::clusterCall(cl, function(){ load_all("~/repos/flowcy/flowcy")}) ## directory that contains the R package.
##' maxres = get_max_lambda_new(ylist, X, numclust, verbose=TRUE,
##'                             nrep = 4,
##'                             ## Function settings
##'                             parallelize = TRUE,
##'                             cl = cl,
##'                             iimax = 6,
##'                             niter = 1000,
##'                             max_lambda_alpha = 10000,
##'                             tol = 1E-3 ## This doesn't need to be so low here.
##'                             )
get_max_lambda_new <- function(ylist, X, numclust,
                               max_lambda_beta = 4000,
                               max_lambda_alpha = 1000,
                               verbose=FALSE,
                               iimax = 16,
                               parallelize = FALSE,
                               cl = NULL,
                               ...){

  ## Get range of regularization parameters.
  ## res0 = covarem_getrange(ylist=ylist, X=X, numclust=numclust, niter=2)

  ## ################################
  ## ## First option: mclapply() ####
  ## ################################
  ## if(is.null(mc.cores)) mc.cores=iimax
  mc.cores = 1
  facs = sapply(1:iimax, function(ii) 2^(-ii)) ## DECREASING order
  print("running the models once")
  if(!parallelize){
    for(ii in 1:iimax){
      if(verbose) printprogress(ii, iimax, fill=TRUE)
      res = covarem(ylist = ylist,
                    X = X,
                    numclust = numclust,
                    pie_lambda = max_lambda_alpha * facs[ii],
                    mean_lambda = max_lambda_beta * facs[ii],
                    verbose=TRUE,
                    ...)

      ## Check zero-ness
      toler = 1E-8
      sum_nonzero_alpha = sum(res$alpha[,-1] > toler)
      sum_nonzero_beta = sum(unlist(lapply(res$beta, function(cf){ sum(cf[-1,] > toler) })))
      if(sum_nonzero_alpha + sum_nonzero_beta != 0 & ii==1){
        stop(paste0("Max lambdas: ", max_lambda_beta, " and ", max_lambda_alpha,
                  " were too small as maximum reg. values. Go up and try again!!"))
      } else {
          myfac = facs[ii-1]
          return(c(max_lambda_beta * myfac, max_lambda_alpha * myfac))
      }
    }
  } else{
    assert_that(!is.null(cl))
    reslist = parallel::parLapplyLB(cl, 1:iimax, function(ii){
      if(verbose) printprogress(ii, iimax, fill=TRUE)
      res = covarem(ylist = ylist,
                    X = X,
                    numclust = numclust,
                    pie_lambda = max_lambda_alpha * facs[ii],
                    mean_lambda = max_lambda_beta * facs[ii],
                    verbose=TRUE,
                    ...)
      return(res[c("alpha", "beta", "mean_lambda", "pie_lambda")])
    })

    ## Then, filter the grid for the number of zeros.
    print("filtering the results")
    toler = 1E-8
    allzero = rep(NA, iimax)
    for(ii in 1:iimax){
      sum_nonzero_alpha = sum(reslist[[ii]]$alpha[,-1] > toler)
      sum_nonzero_beta = sum(unlist(lapply(reslist[[ii]]$beta, function(cf){ sum(cf[-1,] > toler) })))
      allzero[ii] = print(sum_nonzero_alpha + sum_nonzero_beta == 0)
    }
    if(any(allzero)){
      ## Get the SMALLEST coeff such that all zero coefficients were found, and return.
      myfac = facs[max(which(allzero))]
      return(c(max_lambda_beta * myfac, max_lambda_alpha * myfac))
    } else {
      stop(paste0("Max lambdas: ", max_lambda_beta, " and ", max_lambda_alpha,
                  " were too small as maximum reg. values. Go up and try again!!"))
    }
  }

  ## ## #################################################
  ## ## ## New: Do this on a 2d grid; using parLapply. ##
  ## ## #################################################
  ## res0 = list(max_lambda_beta = 2000, max_lambda_alpha = 20)
  ## facs = sapply(1:8, function(ii) 2^(-ii))

  ## ## Do the full 2d thing.
  ## iimat = expand.grid(1:8, 1:8)
  ## iimax = nrow(iimat)
  ## cl1 = get_cl()
  ## reslist = parLapplyLB(cl1, 1:iimax, function(ii){
  ##   fac1 = 2000 * facs[iimat[ii,1]]
  ##   fac2 = 20 * facs[iimat[ii,2]]
  ##   res = covarem(ylist = ylist,
  ##                 X = X,
  ##                 numclust = numclust,
  ##                 mean_lambda = res0$max_lambda_beta,
  ##                 pie_lambda = res0$max_lambda_alpha,
  ##                 ...)
  ##   return(res[c("alpha", "beta", "mean_lambda", "pie_lambda")])
  ## })

  ## ## Then, filter the grid for the number of zeros.
  ## tol = 1E-8
  ## zeromat = matrix(0, nrow = iimax, ncol = 2)
  ## for(ii in 1:iimax){
  ##   sum_zero_alpha = sum(reslist[[ii]]$alpha <= tol)
  ##   sum_zero_beta = sum(unlist(lapply(res$beta, function(cf){ sum(cf[-1,]==0) })) <= tol)
  ##   zeromat[ii,] = c(sum_zero_alpha, sum_zero_beta)
  ## }

  ## allzero = apply(zeromat, 1, function(myrow){
  ##   all(myrow==0)
  ## })

  ## if(any(allzero)){
  ##   ## Randomly choose a guy out of these, and return.
  ##   inds = iimat[sample(which(allzero), 1),]
  ##   return(c(2000 * facs[inds[1]], 20 * facs[inds[2]]))
  ## }
}
