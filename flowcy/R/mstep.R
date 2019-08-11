##' @param resp Is an (T x nt x K) array.
##' @param X Covariate matrix (T x dimdat).
##' @return The multinomial logit model coefficients. A matrix of dimension (K x
##'   (p+1)).
Mstep_alpha <- function(resp, X, numclust, lambda = 0, alpha = 1,
                        refit = FALSE,
                        sel_coef = NULL
                        ){

  TT = nrow(X)
  p = ncol(X)

  ## Calculate the summed responsibilities
  resp.sum = t(sapply(resp, colSums)) ## (T x numclust)

  stopifnot(dim(resp) == c(TT, numclust))

  ## Fit the model, possibly with some regularization
  if(!refit){
    fit = glmnet::glmnet(x = X, y = resp.sum, family = "multinomial",
                         alpha = alpha, lambda = lambda)
    ## The coefficients for the multinomial logit model are (K x (p+1))
    alphahat = do.call(rbind, lapply(coef(fit), t))
    stopifnot(all(dim(alphahat) == c(numclust, (p+1))))

  } else {
    Xa = cbind(1, X)
    alphahat = cvxr_multinom_new(X = Xa, y = resp.sum, lambda = 0,
                                 sel_coef = sel_coef$alpha,
                                 exclude = 1)
    alphahat[which(abs(alphahat) < 1E-8, arr.ind = TRUE)] = 0
    alphahat = t(as.matrix(alphahat))
    stopifnot(all(dim(alphahat) == c(numclust, (p + 1))))
  }


  ## While you're at it, calculate the fitted values (\pi) as well:
  piehatmat = as.matrix(exp(cbind(1, X) %*% t(alphahat)))
  piehat = piehatmat / rowSums(piehatmat)

  stopifnot(all(dim(piehat) == c(TT,numclust)))
  stopifnot(all(piehat >= 0))
  ## This should be dimension (T x K)

  return(list(pie = piehat, alpha = alphahat))##, fit=fit))
}



##' Estimates sigma.
##' @param mn are the fitted means
Mstep_sigma_covar <- function(resp, ylist, mn, numclust){
  TT = length(ylist)
  ntlist = sapply(ylist, nrow)
  dimdat = ncol(ylist[[1]])

  ## Set up empty residual matrix (to be reused)
  resid.rows = matrix(NA, nrow = sum(ntlist), ncol=dimdat)
  cs = c(0, cumsum(ntlist))
  vars = list()
  for(iclust in 1:numclust){

    ## Calculate all the residuals, and weight them.
    for(tt in 1:TT){
      irows = (cs[tt]+1):cs[tt+1]
      row.to.subtract = mn[tt,,iclust]
      resid.rows[irows,] = sqrt(resp[[tt]][,iclust]) *
        sweep(ylist[[tt]], 2, rbind(row.to.subtract))
    }

    ## Calculate the covariance matrix
    resp.grandsum = sum(unlist(sapply(1:TT, function(tt){ resp[[tt]][,iclust] })))
    myvar = crossprod(resid.rows, resid.rows) / resp.grandsum
    vars[[iclust]] = myvar
  }
  rm(resid.rows)

  ## reformat
  sigma = abind::abind(lapply(1:TT, function(tt){
    (abind::abind(vars,along=0))
  }), along=0)
  ## stopifnot(all(dim(sigma) == c(TT, numclust, dimdat, dimdat)))

  ## Return |K| of them.
  return(sigma)
}


## Only works for positive semidefinite matrices that are diagonalizable (no
## normal Jordan forms, etc.)
##' Given a matrix positive definite matrix a, computes a^{-1/2}.
mtsqrt_inv <- function(a){
  a.eig <- eigen(a)
  a.sqrt <- a.eig$vectors %*% diag(1 / sqrt(a.eig$values)) %*% t(a.eig$vectors)
}

##' The M step of beta, using a particular lasso regression formulation.
##' @param maxdev The desired maximum radius of the fitted means from beta0k.
##' @param sel_coef Sparsity pattern. Only active when refit=TRUE.
##' @return Result of M step; a |numclust| length list of (p+1)x(d) matrices,
##'   each containing the estimated coefficients for the mean estimation.
Mstep_beta <- function(resp, ylist, X, mean_lambda=0, sigma, numclust,
                       refit = NULL,
                       sel_coef = NULL,
                       maxdev = NULL
                       ){

  ## Preliminaries
  TT = length(ylist)
  numclust = ncol(resp[[1]])
  dimdat = ncol(ylist[[1]])
  ntlist = sapply(ylist, nrow)
  p = ncol(X)

  Xa = cbind(1, X)

  ## Setup
  manip_obj = manip(ylist, Xa, resp, sigma, numclust)
  Xtildes = manip_obj$Xtildes
  yvecs = manip_obj$yvecs

  ## Intercepts are to be excluded from penalization.
  exclude.from.penalty = (0:(dimdat-1))*(ncol(Xa)) + 1

  ## Obtain fitted beta, separately by cluster.
  results = lapply(1:numclust, function(iclust){

    ## Give the glmnet function pre-calculated Y's and X's.
    if(is.null(maxdev) & !refit){
      ## if(refit) stop("NOT WRITTEN YET") ## Remove when done writing this
      res =  solve_lasso(x = Xtildes[[iclust]], y = yvecs[[iclust]], lambda = mean_lambda,
                         intercept = FALSE,
                         exclude.from.penalty = exclude.from.penalty)
      ## Unravel to obtain the coef and fitted response
      betahat = matrix(res$b, ncol=dimdat)
    } else {
      ## If applicable, restrict the fitted means to be close to \beta0k.
      betahat = cvxr_lasso_newer(X = Xtildes[[iclust]],
                                 Xorig = X,
                                 y = yvecs[[iclust]],
                                 lambda = mean_lambda,
                                 exclude.from.penalty = exclude.from.penalty,
                                 maxdev = maxdev,
                                 dimdat = dimdat,
                                 numclust = numclust,
                                 refit = refit,
                                 sel_coef = sel_coef$beta[[iclust]]
                                 )
      betahat[which(abs(betahat) < 1E-8, arr.ind = TRUE)] = 0
    }
    yhat = Xa %*% betahat
    assert_that(all(dim(betahat) == c(p+1,dimdat)))
    return(list(betahat = betahat, yhat = yhat))
  })

  ## Extract results
  betahats = lapply(results, function(a){a$betahat})
  yhats = lapply(results, function(a){as.matrix(a$yhat)})
  yhats = abind::abind(yhats, along = 0)
  yhats = aperm(yhats, c(2,3,1)) ## This needs to by (T x dimdat x numclust)

  ## Each are lists of length |numclust|.
  return(list(beta = betahats,
              mns = yhats))
}


##' Helper to "manipulate" X and y, to get Xtilde and Ytilde and yvec for a more
##' efficient beta M step (each are |numclust|-length lists, calculated
##' separately for each cluster).
##' @return 3 (or dimdat) |numclust|-length lists.
manip <- function(ylist, X, resp, sigma, numclust){

  ntlist = sapply(ylist, nrow)
  dimdat = ncol(ylist[[1]])
  TT = nrow(X)

  resp.sum = t(sapply(resp, colSums)) ## (T x numclust)
  sigma.inv.halves = array(NA, dim=dim(sigma))
  for(iclust in 1:numclust){
    mymat = mtsqrt_inv(sigma[1,iclust,,])
    for(tt in 1:TT){ sigma.inv.halves[tt,iclust,,] = mymat }
  }

  ## Pre-calculate response and covariates to feed into lasso.
  emptymat = matrix(NA, nrow = dimdat, ncol = TT)
  Ytildes = list()
  for(iclust in 1:numclust){
    cs = c(0, cumsum(ntlist))
    for(tt in 1:TT){
      irows = (cs[tt] + 1):cs[tt + 1]
      wt.y = resp[[tt]][, iclust]  * ylist[[tt]]
      emptymat[,tt] = colSums(wt.y)
    }
    Ytildes[[iclust]] = emptymat
  }
  Xtildes = lapply(1:numclust, function(iclust){
    sigma.inv.halves[1,iclust,,] %x% (sqrt(resp.sum[,iclust]) * X)
  })

  ## Vector
  yvecs = lapply(1:numclust, function(iclust){
    yvec = (1/sqrt(resp.sum[,iclust]) * t(Ytildes[[iclust]])) %*% sigma.inv.halves[1,iclust,,]
    yvec = as.vector(yvec)
  })

  return(list(Xtildes = Xtildes,
              Ytildes = Ytildes,
              yvecs = yvecs))
}
