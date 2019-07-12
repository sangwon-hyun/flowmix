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
    sampled.data = mydata[isamp,,drop=FALSE]
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


##' Initialize  the cluster  centers by  learning clusters  from each  image and
##' aggregating..
##'  @param data  A T-length list of (nt  by 3) datasets.  There should  be T of
##'   such datasets. 3 is actually \code{mulen}.
##' @param numclust Number of clusters (M).
##' @param TT total number of (training) time points.
##' @return An array of dimension (T by numclust by M).
init_mu_warmstart_v3 <- function(data, numclust, TT, tol2, verbose=FALSE){
  muhats = list()
  muhat = init_mu(data, numclust, 1)
  if(verbose){cat(fill=TRUE);cat("Warm start calculation in progress.");  cat(fill=TRUE)}


  ## Collapse /all/ datapoints, and then make some rough centers to start with
  for(tt in TT){
    yt = ylist[[tt]]
    sample(1:nrow(yt), nrow(yt)/TT)
  }
  ## ## Continue here
  ## y.condensed =
  ## y.collapsed = data
  ## data[[1]]

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
  numclust = dim(mu)[2]
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
##' @param warmstart
##' @return List containing the list of mus, pies, and objective values.
driftem <- function(data, mu=NULL, pie=NULL, niter=1000,
                    sigma=NULL, tol1 = 1E-10, tol2 = 1E-4, lam.pi=0, lam.mu=0,
                    s=1E-4, numclust,
                    sigma.fac=1,
                    fix.sigma=FALSE,
                    ## mu.all.times.same=FALSE
                    ## warmstart=FALSE,
                    warmstart = c("none", "seq", "v2"),
                    verbose=FALSE,
                    regularize=FALSE, eval.min=0, eval.max=Inf
                    ){

  ## Define a few things
  if(class(data)!="list") data = list(data)
  dimdat = ncol(data[[1]])
  ntlist = sapply(data, nrow)
  TT = length(data)
  warmstart = match.arg(warmstart)

  ## If missing, generate initial values
  if(is.null(mu)){
    if(warmstart=="v2"){ mu = init_mu_warmstart_v2(data, numclust, TT, tol2, verbose) }
    if(warmstart=="seq"){ mu = init_mu_warmstart(data, numclust, TT, tol2, verbose) }
    if(warmstart=="none"){ mu = init_mu(data, numclust, TT) }
  }
  if(nrow(mu[1,,])!=numclust){
    warning(paste0("Input numclust of ", numclust, "is being replaced by", nrow(mu[1,,])))
    numclust = nrow(mu[1,,])
  }
  if(is.null(pie)) pie = init_pi(numclust, TT)
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

      if(!fix.sigma){
        sigmalist[[iter]] = Mstep_sigma(resp.list, data, numclust, mulist[[iter-1]], TT, dimdat,
                                        regularize=regularize, eval.min=eval.min, eval.max=eval.max)
      } else {
        sigmalist[[iter]] = Mstep_sigma_constant(resp.list, data, numclust, mulist[[iter-1]], TT, dimdat,
                                        regularize=regularize, eval.min=eval.min, eval.max=eval.max)
      }


      ## Recnet change: changed order so that mu is estimated from /updated/ $Sigma$.
      mulist[[iter]] = Mstep_mu_exact(resp.list, numclust, sigmalist[[iter]], mulist[[iter-1]], ## the previous mu is actually not used.
                                      TT, dimdat, lam.mu, data, Xslist=Xslist, DD3.permuted=DD3.permuted,
                                      ntlist=ntlist)

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
              mu_init=mu_init,
              data=data,
              dimdat=dimdat,
              TT=TT,
              numclust=numclust
              ))
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


##' (New  version,  for  now  without  regularization)  Optimize  the  penalized
##' likelihood for pie.  Iterative updates are made to the pie matrix by solving
##' the proximal-type subproblem.
##' @param maxsteps maximum number of inner steps to take.
##' @param pie_init initial value of pie, T by M.
##' @param resp.list TT length list of responsibility matrices, each n by M.
##' @param s step size
##' @param tol numerical tolerance for inner steps
##' @param lam tuning parameter for penalty on pi.
##' @param maxsteps Defaults to 1
Mstep_pi_nopenalty <- function(maxsteps=100, resp.list, TT,
                     pie_init=matrix(NA, ncol=M, nrow=TT), s=1E-1, tol=1E-10,
                     lam){
  resp.collapsed = (do.call(rbind, lapply(resp.list, colSums)))
  return(resp.collapsed/rowSums(resp.collapsed))
}


##' Optimize the penalized likelihood for pie.  Iterative updates are made to the
##' pie matrix by solving the proximal-type subproblem.
##' @param maxsteps maximum number of inner steps to take.
##' @param pie_init initial value of pie, T by M.
##' @param resp.list TT length list of responsibility matrices, each n by M.
##' @param s step size
##' @param tol numerical tolerance for inner steps
##' @param lam tuning parameter for penalty on pi.
##' @param maxsteps Defaults to 1
Mstep_pi <- function(maxsteps=100, resp.list, TT,
                     pie_init=matrix(NA, ncol=M, nrow=TT), s=1E-1, tol=1E-10,
                     lam){
  if(lam==0){
    return(list(pie=Mstep_pi_nopenalty(maxsteps, resp.list, TT, pie_init, s, tol=tol, lam)))
  }

  ## Initialize
  pielist = list()
  objlist = rep(NA, maxsteps)
  istep = 1
  resp.collapsed = (do.call(rbind, lapply(resp.list, colSums)))
  pielist[[istep]] = pie_init
  objlist[istep] = objective_pi(pie_init, resp.collapsed)
  pie_old = pie_init

  ## Iterates
  for(istep in 2:maxsteps){

    ## Update pie
    pie_new = one_pi_update(pie_old, resp.collapsed, TT, lam, s=s)
    pielist[[istep]] = pie_new
    objlist[istep] = objective_pi(pie_new, resp.collapsed)

    ## Stop if converged
    ## print(check_converge(objlist[istep - 1], objlist[istep], tol))
    if(is.nan(objlist[istep])) browser()
    if(check_converge(objlist[istep - 1], objlist[istep], tol)) break
    pie_old = pie_new
  }

  ## Truncate the list
  pielist = pielist[1:istep]
  objlist = objlist[1:istep]
  ## print(pielist[istep])
  return(list(pie=pie_new, pielist=pielist, objlist=objlist))
}

##' Inner update of pie, to find the scaling for u.
##' @param pie pie.
##' @param resp responsibility matrix.
##' @param TT Number of timepoints.
##' @param lam penalty.
##' @param s step size.
one_pi_update <- function(pie, resp, TT, lam, s=1E-5){

  ## Calculate gradient (matrix)
  G = gradient_pi(pie, resp, TT, lam)

  ## Multiply the factor (element-wise)
  pie_new = pie * exp( -s * G - 1)

  ## Scale the rows
  pie_new <- t(apply(pie_new, 1, function(myrow){myrow/sum(myrow)}))
  ## pie_new <- scale(pie_new, center = FALSE,
  ##                  scale = rowSums(pie_new))
  return(pie_new)
}


##' Gradient of objective with respect to pi.
##' @param pie The value of pie at which the gradient is being calculated.
##' @param resp The responsibility matrix.
##' @return Gradient (matrix valued, T by M).
gradient_pi <- function(pie, resp, TT, lam){

  ## Gradient of the likelihood
  grad_loglik = resp/pie

  ## Gradient of the penalty
  ## TODO Figure out the sign of this.
  if(lam > 0 & TT > 1){
    grad_pen = matrix(NA, nrow=nrow(grad_loglik), ncol=ncol(grad_loglik))
    grad_pen[2:(TT-1),] = 2 * (2* pie[2:(TT-1),] - pie[3:TT,] - pie[1:(TT-2),])
    grad_pen[1,] =  2 * (pie[1,] - pie[2,])
    grad_pen[TT,] = 2 * (pie[TT,] - pie[TT-1,])
  } else {
    grad_pen = 0
  }

  ## Entrywise sum of the two gradients
  return(-grad_loglik + lam * grad_pen)
}


##' Gradient of objective with respect to pi.
##' @param  pie   The  value   of  pie   at  which   the  gradient   is  being
##' calculated. Dimension T by M.
##' @param resp The responsibility matrix, dimension T by M
##' @return Gradient (matrix valued, T by M).
objective_pi <- function(pie, resp){
  stopifnot(dim(resp) == dim(pie))
  sum(- resp * log(pie))
}


##' (non-time-varying sigma) Optimize  the penalized Q function  with respect to
##' covariance (closed form solution).
##' @param resp.list List of responsibility matrices.
##' @param  mu array of dimension  T by M by  p. This needs to  be the /updated/
##'   value of mu, in this iteration of EM.
##' @param data T lengthed list of data.
##' @param numclust Number of groups.
##' @param TT total number of timepoints.
##' @param dimdat dimension of data.
##'   @return An  array containing  optimized covariances.   Dimension is  (T by
##'   numclust by dimdat by dimdat).
Mstep_sigma_constant <- function(resp.list, data, numclust, mu, TT, dimdat,
                                 regularize=FALSE, eval.min=0, eval.max=Inf## Added experimentally
                                 ){

  ## Make summed responsibilities (over i)
  resp.sum = t(sapply(resp.list, function(resp){
    apply(resp, 2, sum)
  }))

  ## Double loop across j and t
  sigmalist = array(NA, dim=c(TT, numclust, dimdat, dimdat))
  for(jj in 1:numclust){

    ## Just rbind all the "centered" data across tt!!
    all.centered.dat = lapply(1:TT, function(tt){
      mujt = mu[tt,jj,]
      centered.dat <- sweep(data[[tt]][,], 2, rbind(mujt))
    })
    all.centered.dat = do.call(rbind, all.centered.dat)

    ## Just concatenate all the responsibilities across tt
    all.resp.for.jj = lapply(1:TT, function(tt){ resp.list[[tt]][,jj] })
    all.resp.for.jj = do.call(c, all.resp.for.jj)

    resp.sum.over.tt = sum(resp.sum[,jj])
    fixedsigma = (t(all.centered.dat)  %*%
                            diag(all.resp.for.jj) %*%
                            all.centered.dat) / resp.sum.over.tt

    if(regularize){
      fixedsigma = sigma_regularize(fixedsigma, eval.min, eval.max)
    }

    ## Assign the same fitted covariance to every tt=1:TT.
    for(tt in 1:TT){
      sigmalist[tt,jj,,] = fixedsigma
    }
  }
  return(sigmalist)
}


##' Optimize the  penalized Q function  with respect to covariance  (closed form
##' solution). The covariances are allowed to vary across time.
##' @param resp.list List of responsibility matrices.
##' @param  mu array of dimension  T by M by  p. This needs to  be the /updated/
##'   value of mu, in this iteration of EM.
##' @param data T lengthed list of data.
##' @param numclust Number of groups.
##' @param TT total number of timepoints.
##' @param dimdat dimension of data.
##'   @return An  array containing  optimized covariances.   Dimension is  (T by
##'   numclust by dimdat by dimdat).
Mstep_sigma <- function(resp.list, data, numclust, mu, TT, dimdat,
                        regularize=FALSE, eval.min=0, eval.max=Inf## Added experimentally
                        ){

  ## Make summed responsibilities (over i)
  ## browser()
  ## resp.array = abind::abind(resp.list, along=0)
  ## resp.sum = Reduce("+", lapply(1:n, function(ii) resp.array[,ii,]))
  resp.sum = t(sapply(resp.list, function(resp){
    apply(resp, 2, sum)
  }))

  ## Double loop across j and t
  sigmalist = array(NA, dim=c(TT, numclust, dimdat, dimdat))
  for(jj in 1:numclust){
    for(tt in 1:TT){

      ## Subtract a single row from data, then get outer product
      mujt = mu[tt,jj,]
      centered.dat <- sweep(data[[tt]][,], 2, rbind(mujt))
      sigmahat = (t(centered.dat)  %*%
                            diag(resp.list[[tt]][,jj]) %*%
                  centered.dat) / resp.sum[tt,jj]
      if(regularize){
        sigmahat = sigma_regularize(sigmahat, eval.min, eval.max)
      }
      sigmalist[tt,jj,,] = sigmahat
    }
  }
  return(sigmalist)
}


##' Projects the matrix to have eigenvalue between two values
##' @param eval.min Minimum eigenvalue.
##' @param eval.max Maximum eigenvalue.
sigma_regularize <- function(sigma, eval.min, eval.max){

  ## Temporary setup, for testing
  ## sigma = cov(MASS::mvrnorm(100, mu=c(0,0,0), Sigma=diag(c(2,1,0.5))))
  ## eval.min = 0.8
  ## eval.max = 1.5
  ## End of temporary setup

  ## Obtain eigendecomposition
  decomp = eigen(sigma)
  eigvals = decomp$values

  ## Replace the eigenvalues
  eigvals[eigvals <= eval.min] = eval.min
  eigvals[eigvals >= eval.max] = eval.max

  ## Make the new matrix and return
  newmat = decomp$vectors %*% diag(eigvals) %*% t(decomp$vectors)
  return(newmat)
}



##' Helper  for  makePmat.  It  gives  the location  of  the  x'th  entry  after
##' conversion to the alternative format.
##' @param dimdat dimension of data.
##' @param TT total number of timepoints.
##' @param x index before permutation.
##' @return index /after/ permutation.
myconvert <- function(x, dimdat, TT){
  ((x-1) %% dimdat) * TT + ceiling( x / dimdat )
}

## stopifnot(all.equal(myconvert(1:9,3,3), c(1,4,7,2,5,8,3,6,9)))

##' The  role of  this function  is to  permute columns  so as  to appropriately
##' penalize the successive differences in time (the "alternative" form).
##' @param dimdat dimension of data.
##' @param TT total number of timepoints.
##' @return Permutation matrix.
makePmat <- function(dimdat, TT){
  Pmat = matrix(0,
                nrow=dimdat*TT,
                ncol=dimdat*TT)
  for(tt in 1:(dimdat*TT)){ Pmat[tt, myconvert(tt,dimdat,TT)] = 1 }
  return(Pmat)
}


##' Exactly optimizing the penalized likelihood for mu.   Iterative updates are
##' made to the pie matrix by solving the proximal-type subproblem.
##' @param resp responsibility matrix.
##' @param maxsteps Defaults to 1
##' @param resp.list List of responsibility matrices.
##' @param mu array of dimension T by M by p.
##' @param data T lengthed list of data.
##' @param numclust Number of groups.
##' @param TT total number of timepoints. ##' @param dimdat dimension of data.
##' @param DD3.permuted Only needs to be calculated once.
##' @param Xslist Only needs to be calculated once.
##' @return mu (array, T by M by p). p is \code{dimdat}.
Mstep_mu_exact <- function(resp.list, M, Sigma, mu, TT, dimdat, lam, data,
                           DD3.permuted=NULL,
                           Xslist=NULL,
                           ntlist){



  ## Solve seaparately for each dimension jj=1:M
  muhat.mat = array(0, c(TT, M, dimdat))
  for(jj in 1:M){
    muhat.mat[,jj,] =  Mstep_mu_jj(Sigma, resp.list, jj, ntlist, TT, DD3.permuted, lam, data)
  }
  return(muhat.mat)
}


##' Helper (workhorse) to do a single  Mu step for jj in 1:M, where M  is the number of data
##' dimensions (dimdat).
Mstep_mu_jj <- function(Sigma, resp.list, jj, ntlist, TT, DD3.permuted=NULL, lam, data){


  ## Making a list of Xi vectors; only make them once.
  ## nt = ntlist[[tt]]
  ## Xslist = lapply(1:nt, function(ii){
  ##   Xs = lapply(data[1:TT], function(mydat){ mydat[ii,]})
  ##   Xs = do.call(c, Xs)
  ## })

  ## Form individual inverse covariance
  inv.covariances = inv.covariances.times.X = list()
  for(tt in 1:TT){
    nt = ntlist[[tt]]
    ## Fix a tt
    Sigmahat.ijt = lapply(1:nt, function(ii){
      solve(Sigma[tt,jj,,]) * resp.list[[tt]][ii,jj]})
    inv.covariances[[tt]] = reduce("+", Sigmahat.ijt)
    Xit = lapply(1:nt, function(ii){ as.numeric(data[[tt]][ii,])})
    inv.covariances.times.X[[tt]] = Reduce("+",
                                           mapply(function(a,b){a%*%b},
                                                  Sigmahat.ijt, Xit,
                                                  SIMPLIFY=FALSE))
  }

  rhs = do.call(c, inv.covariances.times.X)
  lhs1 = Matrix::bdiag(inv.covariances)
  if(!is.null(DD3.permuted)){
    lhs2 = lam * nt * DD3.permuted
  } else { lhs2 = 0 }

  ## solve and reformat muhat into (TT by dimdat) matrix
  muhat_jj = solve(lhs1+lhs2, as.numeric(rhs))
  muhat_jj = matrix(muhat_jj, nrow=TT, byrow=TRUE)
  return(muhat_jj)

}


##' Make fused lasso matrix (borrowed from github/genlassoinf).
dual1d_Dmat <- function(m){
  D = matrix(0, nrow = m-1, ncol = m)
  for(ii in 1:(m-1)){
    D[ii,ii] = -1
    D[ii,ii+1] = 1
  }
  return(D)
}

##' (Elementwise, but only  one round, and *without the rest  of the coordinates
##' being updated _for  now_* ; not deleted  yet, just in case it  is useful for
##' coordinatewise  descent).   Optimize  the   penalized  likelihood   for  mu.
##' Iterative updates  are made to the  pie matrix by solving  the proximal-type
##' subproblem.
##' @param resp responsibility matrix.
##' @param maxsteps Defaults to 1
##' @return mu (array, T by M by p). p is \code{dimdat}.
Mstep_mu_elementwise <- function(resp.list, M, Sigma, mu, TT, dimdat, lam, data){

  ## There are no iterations needed, only one round of the closed form solution.
  resp.array = abind::abind(resp.list, along=0)
  resp.sum = Reduce("+", lapply(1:n, function(ii) resp.array[,ii,]))

  opt.beta.mat = array(0, c(TT, M, dimdat))
  ## grad.pen.mat = array(0, c(TT, M, dimdat))
  ## grad.pen.mat[2:(TT-1),,] = 2 * (2 * mu[2:(TT-1),,] - mu[3:TT,,] - mu[(1:(TT-2)),,])
  ## grad.pen.mat[1,,] =  2 * (mu[1,,] - mu[2,,])
  ## grad.pen.mat[TT,,] =  2 * (mu[TT,,] - mu[TT-1,,])

  ## Entrywise solve the QP in closed form.
  for(jj in 1:M){

    ## Experimental outer loop to tt: cycle through until convergence
    converge = FALSE
    tol = 1E-8
    max.iter = 100
    objective = rep(NA, max.iter)
    iter = 1
    objective[iter] = Inf               #Initialize
    opt.beta.mat.list = list()
    opt.beta.mat.list[[1]] = opt.beta

    for(iter in 2:max.iter){

    for(tt in 1:TT){

      ## Experimental
      opt.beta.mat[[iter]] = opt.beta.mat[[iter-1]]
      mu = mu
      ## End of experimental

      inv.sigma.tj = solve(sigma[tt,jj,,])

      ## Calculate various quantities
      Ajt = resp.sum[tt,jj] * inv.sigma.tj
      bjt = Reduce("+",
        lapply(1:n, function(ii){
        -2 * resp.list[[tt]][ii,jj] * (data[[tt]][ii,]) %*% inv.sigma.tj
        })
      )
      Bjt = (if(tt==1 | tt==TT) 1 else 2) * diag(rep(lam, dimdat))
      djt = rep(0, dimdat)
      if (tt>1) djt = djt -2 * (mu[tt-1,jj,])*lam ## Experimental: changing to
      if (tt<TT) djt = djt -2 * (mu[tt+1,jj,])*lam
      opt.beta = solve(2*Ajt + 2*Bjt, t(-bjt - djt))         #This is the optimum
      opt.beta.mat[tt,jj,] = opt.beta

      ## Experimental
      objective_mu <- function(mu){ calculate<-NULL }
      objective[TT*iter + tt] = objective_mu(...)
      ## End of experimental

    }
    converge = (objective[iter*TT]-objective[(iter-1)*TT] < tol)
    if(converge) break
    }

  }
  return(opt.beta.mat)
}
