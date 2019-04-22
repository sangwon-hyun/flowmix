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
init_mu <- function(data, numclust, TT){

  ## Initialize the means
  if(length(data)==1) data = list(data)
  mulist = lapply(data, function(mydata){
    sampled.data = mydata[sample(1:nrow(mydata), numclust),]
    rownames(sampled.data) = paste0("clust", 1:numclust)
    return(sampled.data)
  })
  names(mulist) = 1:TT
  ## if(TT>2){
  muarray = abind::abind(mulist, along=0) ## T by numclust by petal length
  ## } else {
  ##   ## browser()
  ##   muarray = mulist
  ## }
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



##' Update the responsibilities, for time t.
##' @param t time.
##' @param data TT datasets.
##' @param muarray array containing mu.
##' @param sigmaarray array containing sigma.
##' @param pie matrix containing the mixture.
##' @param numclust number of clusters.
##' @return An nt by M matrix of the responsibilities.
update_responsibility <- function(t, data, muarray, sigmaarray, pie, numclust){

  ## ## Calculates conditional probability of the latent membership
  ## cond_prob_for_t <- function(data, t, ii, numclust){


  ##   ## Operating on each dataset
  ##   datapoint = data[[t]][ii,]
  ##   densities = sapply(1:numclust,
  ##                      function(jj){
  ##                        mu = muarray[t,jj,]
  ##                        sigma = sigmaarray[t,jj,,]
  ##                        return(mvtnorm::dmvnorm(datapoint,
  ##                                                mean=mu,
  ##                                                sigma=sigma,
  ##                                                log=FALSE))
  ##                      })
  ##   weighted_densities = densities * pie[t,]
  ##   resp = weighted_densities/sum(weighted_densities)

  ## }
  ## n = nrow(data[[t]])
  ## responsibilities = t(sapply(1:n,
  ##                           function(ii){cond_prob_for_t(data, t, ii,numclust)}))

  ## return(responsibilities)

  ## Get the density of each data point, for each group j=1,..M at a time.
  cond_prob_for_j <- function(data, t, jj, pie){
    mu = muarray[t,jj,]
    sigma = sigmaarray[t,jj,,]
    densities = pie[t,jj] * mvtnorm::dmvnorm(data[[t]],
                                 mean=mu,
                                 sigma=sigma,
                                 log=FALSE)
  }
  responsibilities2 = sapply(1:numclust,
                            function(jj){cond_prob_for_j(data, t, jj, pie)})
  responsibilities2 = 1/rowSums(responsibilities2) * responsibilities2

  return(responsibilities2)
}


##' E step.
Estep <- function(data, TT, mu, sigma, pie, numclust){
  resp.list = lapply(1:TT, function(t){
    resp = update_responsibility(t, data, mu, sigma,
                                 pie, numclust)
    colnames(resp) = paste0("class", 1:numclust)
    return(resp)
  })
  names(resp.list) = paste0("time", 1:TT)
  return(resp.list)
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
  grad_pen = matrix(NA, nrow=nrow(grad_loglik), ncol=ncol(grad_loglik))
  grad_pen[2:(TT-1),] = 2 * (2* pie[2:(TT-1),] - pie[3:TT,] - pie[1:(TT-2),])
  grad_pen[1,] =  2 * (pie[1,] - pie[2,])
  grad_pen[TT,] = 2 * (pie[TT,] - pie[TT-1,])

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


##' Optimize the  penalized Q function  with respect to covariance  (closed form
##' solution).
##' @param resp.list List of responsibility matrices.
##' @param mu array of dimension T by  M by p. (I'm almost sure that) this needs
##'   to be the /updated/ value of mu, in this iteration of EM.
##' @param data T lengthed list of data.
##' @param numclust Number of groups.
##' @param TT total number of timepoints.
##' @param dimdat dimension of data.
##'  @return  An array  containing  optimized  covariances.  Dimension is  T  by
##'   numclust by dimdat by dimdat.
Mstep_sigma <- function(resp.list, data, numclust, mu, TT, dimdat){

  ## Make summed responsibilities (over i)
  ## browser()
  ## resp.array = abind::abind(resp.list, along=0)
  ## resp.sum = Reduce("+", lapply(1:n, function(ii) resp.array[,ii,]))
  resp.sum = t(sapply(resp.list, function(resp){
    apply(resp,2,sum)
  }))

  ## Double loop across j and t
  sigmalist = array(NA, dim=c(TT, numclust, dimdat, dimdat))
  for(jj in 1:numclust){
    for(tt in 1:TT){

      ## Subtract a single row from data, then get outer product
      mujt = mu[tt,jj,]
      centered.dat <- sweep(data[[tt]][,], 2, rbind(mujt))
      centered.dat <- diag(sqrt(resp.list[[tt]][,jj])) %*% centered.dat
      sigmalist[tt,jj,,] = (t(centered.dat) %*% centered.dat) / resp.sum[tt,jj]
    }
  }
  return(sigmalist)
}


##' (wrong version,  delete shortly) Optimize  the penalized likelihood  for mu.
##' Iterative updates  are made to the  pie matrix by solving  the proximal-type
##' subproblem.
##' @param resp responsibility matrix.
##' @param maxsteps Defaults to 1
##' @param resp.list List of responsibility matrices.
##' @param mu array of dimension T by M by p.
##' @param data T lengthed list of data.
##' @param numclust Number of groups.
##' @param TT total number of timepoints.
##' @param dimdat dimension of data.
##' @return mu (array, T by M by p). p is \code{dimdat}.
Mstep_mu <- function(resp.list, M, Sigma, mu, TT, dimdat, lam,data){

  ## There are no iterations needed, only one round of the closed form solution.
  resp.array = abind::abind(resp.list, along=0)
  resp.sum = Reduce("+", lapply(1:n, function(ii) resp.array[,ii,]))

  opt.beta.mat = array(0, c(TT, M, dimdat))

  ## Entrywise solve the QP in closed form.
  for(jj in 1:M){
    for(tt in 1:TT){
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
      if (tt>1) djt = djt -2 * (mu[tt-1,jj,])*lam
      if (tt<TT) djt = djt -2 * (mu[tt+1,jj,])*lam
      opt.beta = solve(2*Ajt + 2*Bjt, t(-bjt - djt))         #This is the optimum
      opt.beta.mat[tt,jj,] = opt.beta
    }
  }
  return(opt.beta.mat)
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
    muhat.mat[,jj,] =  Mstep_mu_jj(Sigma, resp.list, jj, ntlist, TT, DD3.permuted, lam)
  }
  return(muhat.mat)
}


##' Helper (workhorse) to do a single  Mu step for jj in 1:M, where M  is the number of data
##' dimensions (dimdat).
Mstep_mu_jj <- function(Sigma, resp.list, jj, ntlist, TT, DD3.permuted, lam){


  ## Making a list of Xi vectors; only make them once.
  ## nt = ntlist[[tt]]
  ## Xslist = lapply(1:nt, function(ii){
  ##   Xs = lapply(data[1:TT], function(mydat){ mydat[ii,]})
  ##   Xs = do.call(c, Xs)
  ## })

  ## Form individual inverse covariance
  ## inv.covariances = sapply(1:nt, function(ii){

  ## browser()
  inv.covariances = inv.covariances.times.X = list()
  for(tt in 1:TT){
    nt = ntlist[[tt]]
    ## Fix a tt
    Sigmahat.ijt = lapply(1:nt, function(ii){
      solve(Sigma[tt,jj,,]) * resp.list[[tt]][ii,jj]})
    inv.covariances[[tt]] = Reduce("+", Sigmahat.ijt)
    Xit = lapply(1:nt, function(ii){ as.numeric(data[[tt]][ii,])})
    inv.covariances.times.X[[tt]] = Reduce("+",
                                           mapply(function(a,b){a%*%b},
                                                  Sigmahat.ijt, Xit,
                                                  SIMPLIFY=FALSE))

  }
  lhs1 = do.call(Matrix::bdiag, inv.covariances)
  lhs2 = lam * nt * DD3.permuted
  rhs = do.call(c, inv.covariances.times.X)

  lhs1 = Matrix::bdiag(inv.covariances)
  lhs2 = lam * nt * DD3.permuted
  muhat_jj = solve(lhs1+lhs2, as.numeric(rhs))
  ## print(muhat_jj)

  ## reformat muhat
  muhat_jj = matrix(muhat_jj, nrow=TT, byrow=TRUE)
  ## cs = cumsum(ntlist)
  ## browser()
  ## Map(function(a,b){a:b}, c(0,cs[-length(cs)]) + 1, cs[1:(length(cs))])


  return(muhat_jj)

  ## Reformat it back to (TT by dimdat)
  muhat_jj = matrix(muhat_jj, nrow=TT, byrow=TRUE)



    Xit = lapply(1:TT, function(tt){ as.numeric(data[[tt]][ii,])})
    Sigmahat.ijt.times.X = mapply(function(a,b){a%*%b},
                                      Sigmahat.ijt, Xit, SIMPLIFY=FALSE)
    Sigmahat.ijt.times.X.sum = Reduce("+",Sigmahat.ijt.times.X)
    inv.covariances.times.X[[ii]] = Sigmahat.ijt.times.X.sum


## do.call(Matrix::bdiag, inv.covariances)
## do.call(c, inv.covariances.times.X)
##     ## return(Matrix::bdiag(Sigmahat.ijt.sum))

##   ## Old:
##   inv.covariances = sapply(1:nt, function(ii){
##     Sigmahat.ijt = lapply(1:TT, function(tt){
##       solve(Sigma[tt,jj,,])*resp.list[[tt]][ii,jj]})
##     Sigmahat.ijt.sum = Reduce("+", Sigmahat.ijt)
##     return(Matrix::bdiag(Sigmahat.ijt.sum))
##   })


  ## Set up and solve quadratic equation
  ## rhs = Reduce("+", mapply(function(a,b){a%*%b}, inv.covariances, Xslist))
  ## rhs = do.call(c, inv.covariances.times.X)
  lhs1 = Matrix::bdiag(inv.covariances)
  lhs2 = lam * nt * DD3.permuted
  muhat_jj = solve(lhs1+lhs2, as.numeric(rhs))

  ## Reformat it back to (TT by dimdat)
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

##' Optimize the penalized likelihood for mu.  Iterative updates are made to the
##' pie matrix by solving the proximal-type subproblem.
##' @param resp responsibility matrix.
##' @param maxsteps Defaults to 1
##' @return mu (array, T by M by p). p is \code{dimdat}.
Mstep_mu_new <- function(resp.list, M, Sigma, mu, TT, dimdat, lam, data){

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


##' Objectives.
##' @param pie matrix of mixture proportions, T by M.
##' @param mu array of dimension T by M by p.
##' @param data TT lengthed list of data
##' @param sigma array of dimension T by M by p by p.
objective_overall <- function(mu, pie, sigma, data){

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
##' @return List containing the list of mus, pies, and objective values.
driftem <- function(data, mu=NULL, pie=NULL, niter=1000,
                    sigma=NULL, TT, tol1 = 1E-10, tol2 = 1E-4, lam1, lam2, s, numclust){

  ## If missing, generate initial values
  if(is.null(pie)) pie = init_pi(numclust, TT)
  if(is.null(mu)) mu = init_mu(data, numclust, TT)
  if(is.null(sigma)) sigma = init_sigma(data, numclust, TT, fac=1)

  ## Initialize
  dimdat = ncol(data[[1]])
  ntlist = sapply(data, nrow)
  pielist = mulist = sigmalist = list()
  objectives = c()
  pielist[[1]] = pie
  mulist[[1]] = mu
  sigmalist[[1]] = sigma
  objectives[[1]] = objective_overall(mulist[[1]], pielist[[1]], sigmalist[[1]], data)
  if(class(data[[1]])!="matrix"){  data = lapply(data, as.matrix) }

  ## Across jj, the dimension of the data shouldn't change
  ## But across TT, the dimension of the data could very well change
  ## In this case, Xslist doesn't really make sense as written?
  ## The matrix conversion doesn't make sense as well!

  ## Make t(tilde D3) times tilde D3
  D = dual1d_Dmat(TT)
  D3 = Matrix::bdiag(D,D,D)
  D3 = Matrix::bdiag(lapply(1:dimdat, function(mydim){D}))
  D3.permuted =  D3 %*% makePmat(dimdat, TT)
  DD3.permuted = t(D3.permuted) %*% D3.permuted

  for(iter in 2:niter){
    printprogress(iter, niter)

    ## E step, done on each image separately.
    resp.list = Estep(data, TT, mulist[[iter-1]], sigma, pielist[[iter-1]], numclust)

    ## M step: calculate mu / sigma / pi.
    pielist[[iter]] = Mstep_pi(10000, resp.list, TT, pielist[[iter-1]], s, tol=tol1, lam1)$pie   ## Mstep_pi_nopenalty()
    mulist[[iter]] = Mstep_mu_exact(resp.list, numclust, sigmalist[[iter-1]], mulist[[iter-1]],
                          TT, dimdat, lam2, data, Xslist=Xslist, DD3.permuted=DD3.permuted,
                          ntlist=ntlist)
    sigmalist[[iter]] = Mstep_sigma(resp.list, data, numclust, mulist[[iter]], TT, dimdat) # experimental
    objectives[iter] = objective_overall(mulist[[iter]], pielist[[iter]], sigmalist[[iter]], data)

    if(check_converge(objectives[iter - 1],
                      objectives[iter],
                      tol=tol2)) break
  }
  return(list(objectives = objectives[1:iter],
              pielist = pielist[1:iter],
              mulist = mulist[1:iter],
              sigmalist = sigmalist[1:iter],
              final.iter = iter,
              lam1=lam1,
              lam2=lam2))
  ## Consider adding the data object and other configurations in here.
}
