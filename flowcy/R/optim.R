##' Initialize weight matrix
##' @param numclust Number of clusters (M).
##' @param T total number of (training) time points.
##' @return A T by M matrix \code{pi}.
init_pi <- function(numclust, T){
  pie =  matrix(1/numclust,
               nrow=T,
               ncol=numclust)
}

##' Initialize the cluster centers (naively).
##'  @param data  A T-length list of (nt  by 3) datasets.  There should  be T of
##'   such datasets. 3 is actually \code{mulen}.
##' @param numclust Number of clusters (M).
##' @param T total number of (training) time points.
##' @return An array of dimension (T by numclust by M).
init_mu <- function(data, numclust, T){

  ## Initialize the means
  if(length(data)==1) data = list(data)
  mulist = lapply(data, function(mydata){
    sampled.data = mydata[sample(1:nrow(mydata), numclust),]
    rownames(sampled.data) = paste0("clust", 1:numclust)
    return(sampled.data)
  })
  names(mulist) = 1:T
  ## if(T>2){
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
##' @param T total number of (training) time points.
##' @return  An array (T  by M  by dimdat by  dimdat) containing the  (dimdat by
##'   dimdat) covariances.
init_sigma <- function(data, numclust, T, fac=1){

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
##' @param data T datasets.
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
Estep <- function(data, T, mu, sigma, pie, numclust){
  resp.list = lapply(1:T, function(t){
    resp = update_responsibility(t, data, mu, sigma,
                                 pie, numclust)
    colnames(resp) = paste0("class", 1:numclust)
    return(resp)
  })
  names(resp.list) = paste0("time", 1:T)
  return(resp.list)
}


##' (New  version,  for  now  without  regularization)  Optimize  the  penalized
##' likelihood for pie.  Iterative updates are made to the pie matrix by solving
##' the proximal-type subproblem.
##' @param maxsteps maximum number of inner steps to take.
##' @param pie_init initial value of pie, T by M.
##' @param resp.list T length list of responsibility matrices, each n by M.
##' @param s step size
##' @param tol numerical tolerance for inner steps
##' @param lam tuning parameter for penalty on pi.
##' @param maxsteps Defaults to 1
Mstep_pi_new <- function(maxsteps=100, resp.list, T,
                     pie_init=matrix(NA, ncol=M, nrow=T), s=1E-1, tol=1E-10,
                     lam){
  ## Initialize
  pielist = list()
  objlist = rep(NA, maxsteps)
  istep = 1
  resp.collapsed = (do.call(rbind, lapply(resp.list, colSums)))
  pie = resp.collapsed/rowSums(resp.collapsed)
  return(list(pie=pie))
}


##' Optimize the penalized likelihood for pie.  Iterative updates are made to the
##' pie matrix by solving the proximal-type subproblem.
##' @param maxsteps maximum number of inner steps to take.
##' @param pie_init initial value of pie, T by M.
##' @param resp.list T length list of responsibility matrices, each n by M.
##' @param s step size
##' @param tol numerical tolerance for inner steps
##' @param lam tuning parameter for penalty on pi.
##' @param maxsteps Defaults to 1
Mstep_pi <- function(maxsteps=100, resp.list, T,
                     pie_init=matrix(NA, ncol=M, nrow=T), s=1E-1, tol=1E-10,
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
    pie_new = one_pi_update(pie_old, resp.collapsed, T, lam, s=s)
    pielist[[istep]] = pie_new
    objlist[istep] = objective_pi(pie_new, resp.collapsed)

    ## Stop if converged
    ## print(check_converge(objlist[istep - 1], objlist[istep], tol))
    print(istep)
    print(head(objlist))
    if(is.nan(objlist[istep])) browser()
    if(check_converge(objlist[istep - 1], objlist[istep], tol)) break
    pie_old = pie_new
  }

  ## Truncate the list
  pielist = pielist[1:istep]
  objlist = objlist[1:istep]
  print(pielist[[istep]])
  return(list(pie=pie_new, pielist=pielist, objlist=objlist))
}

##' Inner update of pie, to find the scaling for u.
##' @param pie pie.
##' @param resp responsibility matrix.
##' @param T Number of timepoints.
##' @param lam penalty.
##' @param s step size.
one_pi_update <- function(pie, resp, T, lam, s=1E-5){

  ## Calculate gradient (matrix)
  G = gradient_pi(pie, resp, T, lam)

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
gradient_pi <- function(pie, resp, T, lam){

  ## Gradient of the likelihood
  grad_loglik = resp/pie

  ## Gradient of the penalty
  ## TODO Figure out the sign of this.
  grad_pen = matrix(NA, nrow=nrow(grad_loglik), ncol=ncol(grad_loglik))
  grad_pen[2:(T-1),] = 2 * (2* pie[2:(T-1),] - pie[3:T,] - pie[1:(T-2),])
  grad_pen[1,] =  2 * (pie[1,] - pie[2,])
  grad_pen[T,] = 2 * (pie[T,] - pie[T-1,])

  ## Entrywise sum of the two gradients
  return(-grad_loglik + lam * grad_pen)
}


##' Gradient of objective with respect to pi.
##'   @param  pie   The  value   of  pie   at  which   the  gradient   is  being
##'   calculated. Dimension T by M.
##' @param resp The responsibility matrix, dimension T by M
##' @return Gradient (matrix valued, T by M).
objective_pi <- function(pie, resp){
  stopifnot(dim(resp) == dim(pie))
  sum(- resp * log(pie))
}


##' Optimize the  penalized Q function  with respect to covariance  (closed form
##' solution).
##' @param resp.list List of responsibility matrices.
##' @param mu array of dimension T by M by p.
##' @param data T lengthed list of data.
##' @param numclust Number of groups.
##' @param T total number of timepoints.
##' @param dimdat dimension of data.
##' @return An array containing optimized  covariances. Dimension is T by numclust by dimdat
##'   by dimdat.
Mstep_sigma <- function(resp.list, data, numclust, mu, T, dimdat){

  ## Make summed responsibilities (over i)
  resp.array = abind::abind(resp.list, along=0)
  resp.sum = Reduce("+", lapply(1:n, function(ii) resp.array[,ii,]))

  ## Double loop across j and t
  sigmalist = array(NA, dim=c(T, numclust, dimdat, dimdat))
  ## data = lapply(data, as.matrix)
  for(jj in 1:numclust){
    for(tt in 1:T){

      ## Subtract a single row from data, then get outer product
      mujt = mu[tt,jj,]
      centered.dat <- sweep(data[[tt]][,], 2, rbind(mujt))
      centered.dat <- diag(sqrt(resp.list[[tt]][,jj])) %*% centered.dat
      sigmalist[tt,jj,,] = (t(centered.dat) %*% centered.dat) / resp.sum[tt,jj]
    }
  }
  return(sigmalist)
}


##' Optimize the penalized likelihood for mu.  Iterative updates are made to the
##' pie matrix by solving the proximal-type subproblem.
##' @param resp responsibility matrix.
##' @param maxsteps Defaults to 1
##' @return mu (array, T by M by p). p is \code{dimdat}.
Mstep_mu <- function(resp.list, M, Sigma, mu, T, dimdat, n, lam,data){

  ## There are no iterations needed, only one round of the closed form solution.
  resp.array = abind::abind(resp.list, along=0)
  resp.sum = Reduce("+", lapply(1:n, function(ii) resp.array[,ii,]))

  opt.beta.mat = array(0, c(T, M, dimdat))
  ## grad.pen.mat = array(0, c(T, M, dimdat))
  ## grad.pen.mat[2:(T-1),,] = 2 * (2 * mu[2:(T-1),,] - mu[3:T,,] - mu[(1:(T-2)),,])
  ## grad.pen.mat[1,,] =  2 * (mu[1,,] - mu[2,,])
  ## grad.pen.mat[T,,] =  2 * (mu[T,,] - mu[T-1,,])

  # ## Experimental: trying a pure ridge penalty
  ## opt.beta.mat = grad.pen.mat = array(0, c(T, M, dimdat))
  ## grad.pen.mat = 2 *mu
  ## ## End of experimental.

  ## Entrywise solve the QP in closed form.
  for(jj in 1:M){
    for(tt in 1:T){
      inv.sigma.tj = solve(sigma[tt,jj,,])

      ## Calculate various quantities
      Ajt = resp.sum[tt,jj] * inv.sigma.tj
      bjt = Reduce("+",
        lapply(1:n, function(ii){
        -2 * resp.list[[tt]][ii,jj] * (data[[tt]][ii,]) %*% inv.sigma.tj
        })
      )
      Bjt = (if(tt==1 | tt==T) 1 else 2) * diag(rep(lam, dimdat))
      ## djt = (tt>1) * (mu[tt-1,jj,]) + (tt<T) * (mu[tt+1,jj,])
      djt = rep(0, dimdat)
      if (tt>1) djt = djt -2 * (mu[tt-1,jj,])*lam
      if (tt<T) djt = djt -2 * (mu[tt+1,jj,])*lam
      opt.beta = solve(2*Ajt + 2*Bjt, t(-bjt - djt))         #This is the optimum
      opt.beta.mat[tt,jj,] = opt.beta
    }
  }
  return(opt.beta.mat)
}




##' Objectives.
##' @param pie matrix of mixture proportions, T by M.
##' @param mu array of dimension T by M by p.
##' @param data T lengthed list of data
##' @param sigma array of dimension T by M by p by p.
objective_overall <- function(mu, pie, sigma, data){

  loglikelihood <- function(data, t, mu, sigma, pie){
    dat = data[[t]]
    loglik.all.particles = sum(sapply(1:nrow(dat), function(irow){
      ## One particle's log likelihood
      log.lik.particle = sum(sapply(1:numclust, function(iclust){
        mypie = pie[t,iclust]
        mymu = mu[t,iclust,]
        mysigma = sigma[t,iclust,,]
        return(mypie * mvtnorm::dmvnorm(dat[irow,],
                             mean=mymu,
                             sigma=mysigma,
                             log=TRUE))
      }))
    }))
    return(loglik.all.particles)
  }

  ## Calculate the data likelilhood of one time point
  loglikelihoods = sapply(1:T, function(t){
    loglikelihood(data, t, mu, sigma, pie)
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
##' @param T Number of time points.
##' @param tol1 Numerical tolerance for E step.
##' @param tol2 Numerical tolerance for the objective value.
##' @param lam1 Tuning parameter for pi.
##' @param lam2 Tuning parameter for mu.
##' @param s Step size.
##' @param numclust number of clusters.
##' @return List containing the list of mus, pies, and objective values.
driftem <- function(data, mu, pie, niter=1000, sigma, T, tol1 = 1E-10, tol2 = 1E-4, lam1, lam2, s, numclust){

  ## Initialize
  pielist = mulist = sigmalist = list()
  objectives = c()
  pielist[[1]] = pie
  mulist[[1]] = mu
  objectives[[1]] = objective_overall(mulist[[1]], pielist[[1]], sigma, data)
  if(class(data[[1]])!="matrix"){  data = lapply(data, as.matrix) }

  for(iter in 2:niter){
    printprogress(iter, niter)

    ## E step, done on each image separately.
    resp.list = Estep(data, T, mulist[[iter-1]], sigma, pielist[[iter-1]], numclust)

    ## M step: calculate mu and pi, overall.
    obj1 = Mstep_pi_new(10000, resp.list, T, pielist[[iter-1]], s, tol=tol1, lam1)
    obj2 = Mstep_mu(resp.list, numclust, sigmalist[[iter-1]], mulist[[iter-1]], T, dimdat, n, lam2, data)
    obj3 = Mstep_sigma(resp.list, data, numclust, mu, T, dimdat) # experimental

    pielist[[iter]] = obj1$pie
    mulist[[iter]] = obj2
    sigmalist[[iter]] = obj3            # experimental

    objectives[iter] = objective_overall(mulist[[iter]], pielist[[iter]], sigma, data)
    if(check_converge(objectives[iter - 1],
                      objectives[iter],
                      tol=tol2)) break
  }
  return(list(objectives = objectives[1:iter],
              pielist = pielist[1:iter],
              mulist = mulist[1:iter],
              sigmalist = sigmalist[1:iter],
              final.iter = iter))
}
