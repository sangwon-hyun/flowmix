##' Initialize weight matrix
##' @param numclust Number of clusters (M).
##' @param TT total number of (training) time points.
##' @return A T by M matrix \code{pi}.
init_pi <- function(numclust, TT){
  pie =  matrix(1/numclust,
               nrow=TT,
               ncol=numclust)
}

##' Initialize the cluster centers (naively).
##'  @param data  A T-length list of (nt  by 3) datasets.  There should  be T of
##'   such datasets. 3 is actually \code{mulen}.
##' @param numclust Number of clusters (M).
##' @param TT total number of (training) time points.
##' @return An array of dimension (T by numclust by M).
init_mu <- function(data, numclust, TT, all.times.same=FALSE){

  ## If all times are to start with the same data points,
  if(all.times.same){isamp = sample(1:nrow(data[[1]]), numclust)} #Note, this is using the smallest nt.

  ## Initialize the means
  mulist = lapply(1:TT, function(tt){
    mydata = data[[tt]]
    if(!all.times.same){isamp = sample(1:nrow(mydata), numclust)}
    sampled.data = mydata[isamp,]
    rownames(sampled.data) = paste0("clust", 1:numclust)
    return(sampled.data)
  })
  names(mulist) = 1:TT
  muarray = abind::abind(mulist, along=0) ## T by numclust by petal length
  return(muarray)
}


##' Initialize   the  cluster   centers  by   applying  one-image   GMM  results
##' sequentially as starting points for the next image.
##'  @param data  A T-length list of (nt  by 3) datasets.  There should  be T of
##'   such datasets. 3 is actually \code{mulen}.
##' @param numclust Number of clusters (M).
##' @param TT total number of (training) time points.
##' @return An array of dimension (T by numclust by M).
init_mu_warmstart <- function(data, numclust, TT, tol2, verbose=FALSE){
  muhats = list()
  muhat = init_mu(data, numclust, 1)
  if(verbose){cat(fill=TRUE);cat("Warm start calculation in progress.");  cat(fill=TRUE)}
  for(tt in 1:TT){
    if(verbose){cat(fill=TRUE); printprogress(tt, TT, "warm start time points");   cat(fill=TRUE)}
    res = driftem(data[[tt]], numclust=numclust, mu=muhat, tol2=tol2/TT)
    ## potentially more lenient!
    muhat = res$mu[[res$final.iter]]
    muhats[[tt]] = muhat[1,,]
    if(verbose){cat(fill=TRUE)}
  }
  names(muhats) = 1:TT
  muarray = abind::abind(muhats, along=0)
  return(muarray)
}


##' Initialize  the cluster  centers by  learning clusters  from each  image and
##' aggregating..
##'  @param data  A T-length list of (nt  by 3) datasets.  There should  be T of
##'   such datasets. 3 is actually \code{mulen}.
##' @param numclust Number of clusters (M).
##' @param TT total number of (training) time points.
##' @return An array of dimension (T by numclust by M).
init_mu_warmstart_v2 <- function(data, numclust, TT, tol2, verbose=FALSE){
  muhats = list()
  muhat = init_mu(data, numclust, 1)
  if(verbose){cat(fill=TRUE);cat("Warm start calculation in progress.");  cat(fill=TRUE)}

  ktlist = c()
  centers = list()
  for(tt in 1:TT){
    if(verbose){cat(fill=TRUE); printprogress(tt, TT, "warm start time points");   cat(fill=TRUE)}

    ## Obtain best number of clusters in each image
    kt = get_best_kmean_numclust(data[[tt]])
    obj = kmeans(data[[tt]], kt)
    ktlist[tt] = kt
    centers[[tt]] = obj$centers
  }
  ktmax = max(ktlist)
  allcenters = do.call(rbind, centers)
  final.centers = (kmeans(allcenters, ktmax))$centers
  muhats = (lapply(1:TT, function(tt){final.centers}))
  names(muhats) = 1:TT
  muarray = abind::abind(muhats, along=0)
  return(muarray)
}



##' Initialize the covariances (naively).
##' @param data The (nt by 3) datasets. There should be T of them.
##' @param numclust Number of clusters (M).
##' @param TT total number of (training) time points.
##' @return  An array (T  by M  by dimdat by  dimdat) containing the  (dimdat by
##'   dimdat) covariances.
init_sigma <- function(data, numclust, TT, fac=1){

  ndat = nrow(data[[1]])
  pdat = ncol(data[[1]])

  ## data
  sigmalist = lapply(data, function(mydata){
    ## sampled.data = mydata[sample(1:nrow(mydata), numclust),]
    sigmas = lapply(1:numclust, function(iclust){
      onesigma = diag(fac * rep(1, pdat))
      colnames(onesigma) = paste0("datcol", 1:pdat)
      rownames(onesigma) = paste0("datcol", 1:pdat)
      return(onesigma)
    })
    sigmas = abind::abind(sigmas, along=0)
    return(sigmas)
  })
  ## names(mulist) = 1:T
  sigmaarray = abind::abind(sigmalist, along=0) ## T by numclust by petal length
  return(sigmaarray)
}



##' Update the responsibilities, for time tt.
##' @param tt time.
##' @param data TT datasets.
##' @param muarray array containing mu.
##' @param sigmaarray array containing sigma.
##' @param pie matrix containing the mixture.
##' @param numclust number of clusters.
##' @return An nt by M matrix of the responsibilities.
update_responsibility <- function(tt, data, muarray, sigmaarray, pie, numclust){

  ## Get the density of each data point, for each group j=1,..M at a time.
  cond_prob_for_j <- function(data, tt, jj, pie){
    mu = muarray[tt,jj,]
    sigma = sigmaarray[tt,jj,,]
    densities = pie[tt,jj] * mvtnorm::dmvnorm(data[[tt]],
                                              mean=mu,
                                              sigma=sigma,
                                              log=FALSE)
  }
  responsibilities2 = sapply(1:numclust,
                            function(jj){cond_prob_for_j(data, tt, jj, pie)})
  responsibilities2 = 1/rowSums(responsibilities2) * responsibilities2

  return(responsibilities2)
}


##' E step.
Estep <- function(data, TT, mu, sigma, pie, numclust){
  resp.list = lapply(1:TT, function(tt){
    resp = update_responsibility(tt, data, mu, sigma,
                                 pie, numclust)
    colnames(resp) = paste0("class", 1:numclust)
    return(resp)
  })
  names(resp.list) = paste0("time", 1:TT)
  return(resp.list)
}



##' Objectives.
##' @param pie matrix of mixture proportions, T by M.
##' @param mu array of dimension T by M by p.
##' @param data TT lengthed list of data
##' @param sigma array of dimension T by M by p by p.
objective_overall <- function(mu, pie, sigma, data){
  TT = length(data)
  loglikelihood <- function(data, tt, mu, sigma, pie){
    dat = data[[tt]]
    loglik.all.particles = sum(sapply(1:nrow(dat), function(irow){
      ## One particle's log likelihood
      log.lik.particle = sum(sapply(1:numclust, function(iclust){
        mypie = pie[tt,iclust]
        mymu = mu[tt,iclust,]
        mysigma = sigma[tt,iclust,,]
        return(mypie * mvtnorm::dmvnorm(dat[irow,],
                             mean=mymu,
                             sigma=mysigma,
                             log=TRUE))
      }))
    }))
    return(loglik.all.particles)
  }

  ## Calculate the data likelilhood of one time point
  loglikelihoods = sapply(1:TT, function(tt){
    loglikelihood(data, tt, mu, sigma, pie)
  })

  ## loglikelihood(datalist[[1]], mulist, pielist, sigmalist)
  sum(loglikelihoods)
}


##' The EM to fit the drifting-clusters.
##' @param data Data
##' @param mu Initial means.
##' @param pie Initial pie.
##' @param niter Number of iterations.
##' @param sigma Array containing covariances (initial value).
##' @param TT Number of time points.
##' @param tol1 Numerical tolerance for E step.
##' @param tol2 Numerical tolerance for the objective value.
##' @param lam1 Tuning parameter for pi.
##' @param lam2 Tuning parameter for mu.
##' @param s Step size.
##' @param numclust number of clusters.
##' @param sigma.fac Size of marginal variance in the initial sigmas.
##' Wparam warmstart
##' @return List containing the list of mus, pies, and objective values.
driftem <- function(data, mu=NULL, pie=NULL, niter=1000,
                    sigma=NULL, tol1 = 1E-10, tol2 = 1E-4, lam.pi=0, lam.mu=0,
                    s=1E-4, numclust,
                    sigma.fac=1,
                    fix.sigma=FALSE,
                    ## mu.all.times.same=FALSE
                    ## warmstart=FALSE,
                    warmstart = c("none", "seq", "v2"),
                    verbose=FALSE
                    ){

  ## Define a few things
  if(class(data)!="list") data = list(data)
  dimdat = ncol(data[[1]])
  ntlist = sapply(data, nrow)
  TT = length(data)
  warmstart = match.arg(warmstart)

  ## If missing, generate initial values
  if(is.null(pie)) pie = init_pi(numclust, TT)
  if(is.null(mu)){
    if(warmstart=="v2"){ mu = init_mu_warmstart_v2(data, numclust, TT, tol2, verbose) }
    if(warmstart=="seq"){ mu = init_mu_warmstart(data, numclust, TT, tol2, verbose) }
    if(warmstart=="none"){ mu = init_mu(data, numclust, TT) }
  }
  if(is.null(sigma)) sigma = init_sigma(data, numclust, TT, fac=sigma.fac)

  ## Initialize
  pielist = mulist = sigmalist = list()
  objectives = rep(NA, niter)
  pielist[[1]] = pie
  mulist[[1]] = mu_init = mu
  sigmalist[[1]] = sigma
  objectives[[1]] = objective_overall(mulist[[1]], pielist[[1]], sigmalist[[1]], data)
  if(class(data[[1]])!="matrix"){  data = lapply(data, as.matrix) }
  sound = TRUE

  ## Make t(tilde D3) times tilde D3
  if(TT > 1){
    DD3.permuted = make_DD3_permuted(TT, dimdat)
  } else {
    if(lam.mu!=0 | lam.pi !=0) stop("Can't use regularization if there isonly one time point!")
    DD3.permuted = NULL
  }

  tryCatch({
    for(iter in 2:niter){
      if(verbose)printprogress(iter, niter, "iteration")

      ## E step, done on each image separately.
      resp.list = Estep(data, TT, mulist[[iter-1]], sigmalist[[iter-1]], pielist[[iter-1]], numclust)

      ## M step: calculate mu / sigma / pi.
      pielist[[iter]] = Mstep_pi(10000, resp.list, TT, pielist[[iter-1]], s, tol=tol1, lam.pi)$pie
      mulist[[iter]] = Mstep_mu_exact(resp.list, numclust, sigmalist[[iter-1]], mulist[[iter-1]],
                                      TT, dimdat, lam.mu, data, Xslist=Xslist, DD3.permuted=DD3.permuted,
                                      ntlist=ntlist)

      if(!fix.sigma){
        sigmalist[[iter]] = Mstep_sigma(resp.list, data, numclust, mulist[[iter]], TT, dimdat)
      } else {
        sigmalist[[iter]] = Mstep_sigma_constant(resp.list, data, numclust, mulist[[iter]], TT, dimdat)
      }
      objectives[iter] = objective_overall(mulist[[iter]], pielist[[iter]], sigmalist[[iter]], data)

      if(check_converge_rel(objectives[iter - 1],
                            objectives[iter],
                            tol=tol2)) break
    }
  },  error = function(err) {
    err$message = paste(err$message,"\n(EM computation has been terminated;",
                        " partial results are being returned.)",sep="")
    warning(err)
  })

  ## As a hackish check, truncate further
  ## if(is.null(sigmalist[[iter]]) | is.null(mulist[[iter]]))  iter = iter-1
  if(length(sigmalist)!=iter | length(sigmalist)!=iter )  iter = iter-1 #Alternatively,
                                                                        #final.iter
                                                                        #can be
                                                                        #adjusted.

  ## Trim and return
  return(list(objectives = objectives[1:iter],
              pielist = pielist[1:iter],
              mulist = mulist[1:iter],
              sigmalist = sigmalist[1:iter],
              final.iter = iter,
              lam.pi=lam.pi,
              lam.mu=lam.mu,
              sound=sound,
              mu_init=mu_init))
  ## Consider adding the data object and other configurations in here.
}


##' Wrapper to run EM with progressively small step size for the exponentiated gradient descent.
driftem_wrapper <- function(data, mu=NULL, pie=NULL, niter=1000,
                    sigma=NULL, tol1 = 1E-10, tol2 = 1E-4, lam.pi=0, lam.mu=0,
                    s=1E-4, numclust,
                    sigma.fac=1,
                    fix.sigma=FALSE,
                    mu.all.times.same=FALSE){

  ## Run while sound
  s.curr = s
  sound = FALSE
  while(!sound){
    cat(fill=TRUE)
    res = driftem(data, mu, pie, niter,
                    sigma, tol1, tol2, lam.pi, lam.mu,
                    s.curr, numclust,
                    sigma.fac,
                    fix.sigma,
                    mu.all.times.same)
    print(paste0("running with s=", s.curr))
    cat(fill=TRUE)
    s.curr = s.curr * 0.1
  }
  return(res)
}



##' Making the permuted version of \tilde D \tilde D.
make_DD3_permuted <- function(TT, dimdat){
  D = dual1d_Dmat(TT)
  D3 = Matrix::bdiag(lapply(1:dimdat, function(mydim){D}))
  D3.permuted =  D3 %*% makePmat(dimdat, TT)
  DD3.permuted = t(D3.permuted) %*% D3.permuted
  return(DD3.permuted)
}
