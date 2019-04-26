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
Mstep_sigma_constant <- function(resp.list, data, numclust, mu, TT, dimdat){

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

    ## Assign the same fitted covariance to every tt=1:TT.
    for(tt in 1:TT){
      sigmalist[tt,jj,,] = fixedsigma
    }
  }
  return(sigmalist)
}


##' Optimize the  penalized Q function  with respect to covariance  (closed form
##' solution).
##' @param resp.list List of responsibility matrices.
##' @param  mu array of dimension  T by M by  p. This needs to  be the /updated/
##'   value of mu, in this iteration of EM.
##' @param data T lengthed list of data.
##' @param numclust Number of groups.
##' @param TT total number of timepoints.
##' @param dimdat dimension of data.
##'   @return An  array containing  optimized covariances.   Dimension is  (T by
##'   numclust by dimdat by dimdat).
Mstep_sigma <- function(resp.list, data, numclust, mu, TT, dimdat){

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
      sigmalist[tt,jj,,] = (t(centered.dat)  %*%
                            diag(resp.list[[tt]][,jj]) %*%
                            centered.dat) / resp.sum[tt,jj]
    }
  }
  return(sigmalist)
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
Mstep_mu_jj <- function(Sigma, resp.list, jj, ntlist, TT, DD3.permuted=NULL, lam){


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
    inv.covariances[[tt]] = Reduce("+", Sigmahat.ijt)
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
