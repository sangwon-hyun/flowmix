##' Main function for covariate EM.
##' @param ylist T-length list each containing response matrices of size (nt x 3), which
##'   contains coordinates of  the 3-variate particles, organized  over time (T)
##'   and with (nt) particles at every time.
##' @param X Matrix of size (T x p+1)
##' @param beta.list Each element is a (T x p+1 x 3 x K) array
##' @param alpha.list Each element is a T by p+1 array
##' @param resp.list  Each element is the posterior probabilities  of the latent
##' @param pie.list (T by K)
##'   variable $Z$ given the paramter estimates. Each  element is a T by nt by k
##'   lengthed vectors
covarem <- function(ylist, X, numclust, niter=100, mn=NULL){


  ## some other things to be defined
  ntlist = sapply(ylist, nrow)
  dimdat = ncol(ylist[[1]])
  TT = length(ylist)
  p = ncol(X)

  ## Data structures
  beta = init_beta(TT, p, dimdat, numclust)
  alpha = init_alpha(TT, p)
  ## mn = calc_mu(beta, X, dimdat, numclust) ## (T x dimdat x numclust)
  if(is.null(mn)){
    mn = aperm(init_mu(lapply(ylist, cbind), numclust, TT), c(1,3,2))
    print('here')
  }

  pie = calc_pie(TT, numclust) ## Let's just say it is all 1/K for now.
  sigma = init_sigma(ylist, numclust, TT) ## (T x numclust x dimdat x dimdat)

  ## Initialize alpha and beta
  beta.list = alpha.list = resp.list = sigma.list = pie.list = mn.list = list()
  objectives = rep(NA, niter)
  objectives[1] = -1E20
  beta.list[[1]] = beta
  alpha.list[[1]] = alpha
  mn.list[[1]] = mn
  sigma.list[[1]] = sigma
  pie.list[[1]] = pie
  ## objectives[1] = objectives(...)       #Insert initial values

  for(iter in 2:niter){
    printprogress(iter,niter, "EM iterations.")

    ## Conduct E step
    resp.list[[iter]] <- Estep_covar(mn.list[[iter-1]],
                                     sigma.list[[iter-1]],
                                     pie.list[[iter-1]],
                                     ylist,
                                     numclust,
                                     ntlist,
                                     iter=iter)  ## This should be (T x numclust x dimdat x dimdat)

    ## Conduct M step

    ## 1. Sigma
    sigma.list[[iter]] <- Mstep_sigma_covar(resp.list[[iter]],
                                            ylist,
                                            mn.list[[iter-1]],
                                            numclust)

    ## 2. Alpha
    res.alpha = Mstep_alpha(resp.list[[iter]],
                            X, numclust, iter=iter)
    alpha.list[[iter]] = res.alpha$alpha
    pie.list[[iter]] = res.alpha$pie

    ## 3. Beta
    res.beta = Mstep_beta(resp.list[[iter]], ylist, X)
    beta.list[[iter]] = res.beta$beta
    mn.list[[iter]]    = res.beta$mns

    ## Calculate the objectives
    objectives[iter] = objective_overall_cov(aperm(mn.list[[iter]], c(1,3,2)),
                                             pie.list[[iter]],
                                             sigma.list[[iter]], ylist)

    ## Check convergence
    if(check_converge_rel(objectives[iter-1],
                          objectives[iter], tol=1E-6))break
  }

  return(list(alpha.list=alpha.list[1:iter],
              beta.list=beta.list[1:iter],
              resp.list=resp.list[1:iter],
              mn.list=mn.list[1:iter],
              sigma.list=sigma.list[1:iter],
              pie.list=pie.list[1:iter],
              objectives=objectives[1:iter],
              final.iter=iter))
}


##' Objectives.
##' @param pie matrix of mixture proportions, T by M.
##' @param mu array of dimension T by M by p.
##' @param data TT lengthed list of data
##' @param sigma array of dimension T by M by p by p.
objective_overall_cov <- function(mu, pie, sigma, data){
  TT = length(data)
  loglikelihood_tt <- function(data, tt, mu, sigma, pie){
    dat = data[[tt]]
    ## One particle's log likelihood
    log.lik.allclust = sapply(1:numclust, function(iclust){

    ## loglik.all.particles = sum(sapply(1:nrow(dat), function(irow){
      mydat = dat
      mypie = pie[tt,iclust]
      mymu = mu[tt,iclust,]
      mysigma = matrix(sigma[tt,iclust,,])
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

  -sum(loglikelihoods)
}
