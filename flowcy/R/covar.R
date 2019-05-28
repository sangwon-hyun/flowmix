##' Main function for covariate EM.
##' @param ylist  T-length list each containing response matrices  of size (nt x
##'   3), which contains coordinates of  the 3-variate particles, organized over
##'   time (T) and with (nt) particles at every time.
##' @param X Matrix of size (T x p+1)
##' @param pie.list (T by K)
##' @param mean_lambda lambda for lasso for the mean.
##' @param pie_lambda lambda for lasso for pie.
##' @return  List containing  fitted parameters and  means and  mixture weights,
##'   across algorithm iterations.
covarem <- function(ylist, X, numclust, niter=100, mn=NULL, pie_lambda=0,
                    mean_lambda=0, verbose=FALSE,
                    warmstart = c("none", "rough"), sigma.fac=1){

  ## Setup.
  ntlist = sapply(ylist, nrow)
  dimdat = ncol(ylist[[1]])
  TT = length(ylist)
  p = ncol(X)
  warmstart = match.arg(warmstart)

  ## Initialize.
  beta = init_beta(TT, p, dimdat, numclust)
  alpha = init_alpha(TT, p)
  if(is.null(mn)){
    if(warmstart=="rough"){
      mn = warmstart_covar(ylist, numclust)
    } else if (warmstart=="none"){
      mn = aperm(init_mu(lapply(ylist, cbind), numclust, TT), c(1,3,2))
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
  for(iter in 2:niter){
    if(verbose) printprogress(iter, niter, "EM iterations.", start.time=start.time)

    ## Conduct E step
    resp <- Estep_covar(mn.list[[iter-1]],
                        sigma.list[[iter-1]],
                        pie.list[[iter-1]],
                        ylist,
                        numclust,
                        ntlist)  ## This should be (T x numclust x dimdat x dimdat)

    ## Conduct M step
    ## 1. Alpha
    res.alpha = Mstep_alpha(resp,
                            X, numclust, iter=iter,
                            lambda=pie_lambda)
    alpha.list[[iter]] = res.alpha$alpha
    pie.list[[iter]] = res.alpha$pie

    ## 2. Beta
    res.beta = Mstep_beta_faster_lasso(resp, ylist, X,
                                       mean_lambda=mean_lambda,
                                       sigma.list[[iter-1]])
    beta.list[[iter]] = res.beta$beta
    mn.list[[iter]]    = res.beta$mns

    ## 3. Sigma
    sigma.list[[iter]] <- Mstep_sigma_covar(resp,
                                            ylist,
                                            mn.list[[iter]],
                                            numclust)

    ## Calculate the objectives
    objectives[iter] = objective_overall_cov(aperm(mn.list[[iter]], c(1,3,2)),
                                             pie.list[[iter]],
                                             sigma.list[[iter]],
                                             ylist,
                                             pie_lambda=pie_lambda,
                                             mean_lambda=mean_lambda,
                                             alpha=res.alpha$alpha,
                                             beta=res.beta$beta)

    ## Check convergence
    if(check_converge_rel(objectives[iter-1],
                          objectives[iter], tol=1E-6)) break
  }


  return(list(alpha.list=alpha.list[1:iter],
              alpha.fit=res.alpha$fit,
              beta.list=beta.list[1:iter],
              mn.list=mn.list[1:iter],
              sigma.list=sigma.list[1:iter],
              pie.list=pie.list[1:iter],
              objectives=objectives[1:iter],
              final.iter=iter,
              ## Above is output, below are data/algorithm settings.
              ntlist=ntlist,
              dimdat=dimdat,
              TT=TT,
              p=p,
              numclust=numclust,
              ## ylist=ylist,
              X=X,
              pie_lambda=pie_lambda,
              mean_lambda=mean_lambda
              ))

}


##' Objectives.
##' @param pie matrix of mixture proportions, T by M.
##' @param mu array of dimension T by M by p.
##' @param data TT lengthed list of data
##' @param sigma array of dimension T by M by p by p.
##' @param alpha linear coefficients for regression on (log ratio of) pie.
##' @param beta linear coefficients for regression on mean.
objective_overall_cov <- function(mu, pie, sigma, data, pie_lambda=0, mean_lambda=0, alpha=0,
                                  beta=0){
  TT = length(data)
  numclust = dim(mu)[2] ## Temporary; there must be a better solution for this.
  loglikelihood_tt <- function(data, tt, mu, sigma, pie){
    dat = data[[tt]]
    ## One particle's log likelihood
    log.lik.allclust = sapply(1:numclust, function(iclust){
      mydat = dat
      mypie = pie[tt,iclust]
      mymu = mu[tt,iclust,]
      mysigma = as.matrix(sigma[tt,iclust,,])
      return(mypie * mvtnorm::dmvnorm(mydat,
                                      mean=mymu,
                                      sigma=mysigma,
                                      log=TRUE))
      })
    return(log.lik.allclust)
  }

  ## Calculate the data likelilhood of one time point
  loglikelihoods = sapply(1:TT, function(tt){
    loglikelihood_tt(data, tt, mu, sigma, pie)
  })

  ## TODO: Form and add penalty term.
  l1norm <- function(coef){ sum(abs(coef)) }
  -sum(unlist(loglikelihoods)) + pie_lambda * l1norm(as.numeric(unlist(alpha))) + mean_lambda * l1norm(as.numeric(unlist(beta)))
}



##' Warmstarts for covariate EM.
##' @param ylist list of data.
##' @param numclust number of clusters desired.
warmstart_covar <- function(ylist, numclust){

  dimdat = ncol(ylist[[1]])

  ## Collapse all the data
  all.y = do.call(rbind, ylist)

  ## Do a cheap k-means
  obj = kmeans(all.y, numclust)
  centres = array(NA, dim=c(TT, dimdat, numclust))
  for(tt in 1:TT){
    centres[tt,,] = t(obj$centers)
  }

  ## Repeat it TT times and return it
  return(centres)
}


##' Prediction: given  new X's,  generate a  set of means  and pies  (and return
##' default Sigma)
##' @param res object returned from covariate EM covarem().
predict.covarem <- function(res, newx=NULL){

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

  ## Predict the means (manually).
  beta = res$beta[[res$final.iter]]
  newmn = lapply(1:numclust, function(iclust){
    newx.a  %*%  beta[[iclust]]
  })
  newmn = abind::abind(newmn, along=0)
  newmn = aperm(newmn, c(2,3,1)) ## This needs to by (T x dimdat x numclust)

  ## Predict the pies.
  newpie = predict(res$alpha.fit, newx=newx)[,,1]

  ## Return all three things
  return(list(newmn=newmn,
              newpie=newpie,
              sigma=res$sigma[[res$final.iter]]))
}
