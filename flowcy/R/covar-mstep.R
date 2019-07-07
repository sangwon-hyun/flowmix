
##' @param resp Is an (T x nt x K) array.
##' @param X Covariate matrix (T x dimdat).
##' @return The multinomial logit model coefficients. A matrix of dimension (K x
##'   (p+1)).
Mstep_alpha <- function(resp, X, numclust, lambda=0, alpha=1,
                        refit=FALSE,
                        sel_coef=NULL
                        ## coef_limit=NULL ## Experimental feature
                        ){

  TT = nrow(X)
  p = ncol(X)

  ## Calculate the summed responsibilities
  resp.sum = t(sapply(resp, colSums)) ## (T x numclust)

  stopifnot(dim(resp)==c(TT, numclust))

  ## Fit the model, possibly with some regularization
  if(!refit){
    fit = glmnet::glmnet(x=X, y=resp.sum, family="multinomial",
                         alpha=alpha, lambda=lambda)
    ## The coefficients for the multinomial logit model are (K x (p+1))
    alphahat = do.call(rbind, lapply(coef(fit), t))
    stopifnot(all(dim(alphahat)==c(numclust, (p+1))))

  } else {
    ## Temporary, for the sel_coef feature
    Xa = cbind(1, X)
    alphahat = cvxr_multinom_new(X=Xa, y=resp.sum, lambda=lambda,
                                 sel.coef=sel_coef$alpha,
                                 exclude=1)
    alphahat[which(abs(alphahat) < 1E-8, arr.ind=TRUE)] = 0
    alphahat = t(as.matrix(alphahat))
    stopifnot(all(dim(alphahat)==c(numclust, (p+1))))
    ## End of temporary
  }

  ## While you're at it, calculate the fitted values (\pi) as well:
  piehatmat = as.matrix(exp(cbind(1,X) %*% t(alphahat)))
  piehat = piehatmat / rowSums(piehatmat)

  stopifnot(all(dim(piehat)==c(TT,numclust)))
  stopifnot(all(piehat>=0))
  ## This should be dimension (T x K)

  return(list(pie=piehat, alpha=alphahat))##, fit=fit))
}


##' Essentially a regression.
Mstep_beta <- function(resp, ylist, X){

  Xa = cbind(rep(1,nrow(X)), X)
  numclust = ncol(resp[[1]])

  resp.sum = t(sapply(resp, colSums)) ## (T x numclust)

  ## Gather Xtilde
  wt.X.list = lapply(1:numclust, function(iclust){
    Xk = sqrt(resp.sum[,iclust]) * Xa ## The left vector should operate on the rows
  }) ## K lengthed list of (T x p+1) matrices

  ## Gather Y
  wt.Y.list = lapply(1:numclust, function(iclust){
    Yk.list = lapply(1:TT, function(tt){
      apply(resp[[tt]][,iclust] * ylist[[tt]], 2, sum)
    })
    Yk = do.call(rbind, Yk.list)
  })

  ## Calculate coefficients
  betahats = lapply(1:numclust, function(iclust){
    Xtilde = wt.X.list[[iclust]]
    solve(t(Xtilde) %*% Xtilde,
          t(Xa) %*% wt.Y.list[[iclust]])
  })
  ## |numclust| length list containing matrices each of dimension (p+1 x dimdat)
  ## Possibly think about regularization here.

  yhats = lapply(betahats, function(betahat){
    Xa %*%  betahat ## Double-check this.
  })
  yhats = abind::abind(yhats, along=0)
  yhats = aperm(yhats, c(2,3,1))
  ## |numclust| length list containing matrices each of dimension (T x dimdat)

  ## Each are lists of length |numclust|.
  return(list(beta=betahats,
              mns=yhats))

  ## to do: Maybe manually do a QR decomp instead of solve().
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

    ## ## Can we do it in one swipe? (will it help?)
    ## subtract.row = rbind(1:3)
    ## tt=1
    ## for(tt in 1:TT){
    ##   irows = (cs[tt]+1):cs[tt+1]
    ##   resid.rows[irows,] = matrix(rep(subtract.row, ntlist[tt]), ncol=dimdat, byrow=TRUE)
    ## }
    ## resp[[tt]][,iclust]
    ## ## End of fix (under construction)

    ## ## Calculate all the residuals, and weight them.
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




##' Making it into a row-stacked setup, so that lasso can be applied (eventually
##' combine it back into Mstep_beta as an option..
Mstep_beta_lasso <- function(resp, ylist, X, mean_lambda=0){

  TT = length(ylist)
  numclust = ncol(resp[[1]])
  ntlist = sapply(ylist, nrow)
  Xa = cbind(rep(1,nrow(X)), X)

  ## For each cluster, separately
  results = lapply(1:numclust, function(iclust){
    ## Stack and reweight covariates into regression form.
    Xaa = do.call(rbind, lapply(1:TT, function(tt){
      myresp = resp[[tt]][, iclust]
      Xrep = matrix(rep(Xa[tt,],each=ntlist[tt]),nrow=ntlist[tt])
      (sqrt(myresp)) * Xrep ## Make sure this applies to the /rows/.
    }))

    ## Now, stack y into a long row
    wt.ylist = lapply(1:TT, function(tt){
      myresp = resp[[tt]][, iclust]
      (sqrt(myresp)) * ylist[[tt]]
    })
    yaa = do.call(rbind, wt.ylist)
    assert_that(nrow(yaa) == sum(ntlist))

    ## Then, do a regular regression
    fit = glmnet::glmnet(y=yaa, x=Xaa, lambda=mean_lambda,
                         alpha=1, intercept=FALSE, family="mgaussian")

    betahat = do.call(cbind, coef(fit))[-1,]
    yhat = Xa %*%  betahat ## Double-check this.
    return(list(betahat=betahat, yhat=yhat))
  })
  ## ## I could have also achieved using the weighted lasso with responsibilities
  ## ## as weights. Check this.

  ## Extract restuls
  betahats = lapply(results, function(a){a$betahat})
  yhats = lapply(results, function(a){as.matrix(a$yhat)})
  yhats = abind::abind(yhats, along=0)
  yhats = aperm(yhats, c(2,3,1)) ## This needs to by (T x dimdat x numclust)

  ## Each are lists of length |numclust|.
  return(list(beta=betahats,
              mns=yhats))
}

## Only works for positive semidefinite matrices that are diagonalizable (no
## normal Jordan forms, etc.)
##' Given a matrix positive definite matrix a, computes a^{-1/2}.
mtsqrt_inv <- function(a){
  a.eig <- eigen(a)
  a.sqrt <- a.eig$vectors %*% diag(1 / sqrt(a.eig$values)) %*% t(a.eig$vectors)
}

##' The M step of beta using a faster lasso regression formulation.
Mstep_beta_faster_lasso <- function(resp, ylist, X, mean_lambda=0, sigma, numclust,
                                    ## coef_limit=NULL## Experimental feature
                                    refit=NULL,
                                    sel_coef=NULL
                                    ){

  ## Preliminaries
  TT = length(ylist)
  numclust = ncol(resp[[1]])
  dimdat = ncol(ylist[[1]])
  ntlist = sapply(ylist, nrow)

  Xa = cbind(rep(1, TT), X)

  ## Setup
  resp.sum = t(sapply(resp, colSums)) ## (T x numclust)
  sigma.inv.halves = array(NA, dim=dim(sigma))
  for(iclust in 1:numclust){
    mymat = mtsqrt_inv(sigma[1,iclust,,])
    for(tt in 1:TT){ sigma.inv.halves[tt,iclust,,] = mymat }
  }

  ## Pre-calculate response and covariates to feed into lasso.
  emptymat = matrix(NA, nrow=dimdat, ncol=TT)
  Ytilde = list()
  for(iclust in 1:numclust){
    cs = c(0, cumsum(ntlist))
    for(tt in 1:TT){
      irows = (cs[tt]+1):cs[tt+1]
      wt.y = resp[[tt]][, iclust]  * ylist[[tt]]
      emptymat[,tt] = colSums(wt.y)
    }
    Ytilde[[iclust]] = emptymat
  }
  Xtilde = lapply(1:numclust, function(iclust){
    sigma.inv.halves[1,iclust,,] %x% (sqrt(resp.sum[,iclust]) * Xa)
  })

  ## Intercepts are to be excluded.
  ## penalty.factor = rep(1, ncol(Xa) * dimdat)
  ## penalty.factor = penalty.factor[(0:(dimdat-1))*(ncol(Xa)) + 1]## = 0
  exclude.from.penalty = (0:(dimdat-1))*(ncol(Xa)) + 1## = 0

  results = lapply(1:numclust, function(iclust){

    ## Form y vector
    yvec = (1/sqrt(resp.sum[,iclust]) * t(Ytilde[[iclust]])) %*% sigma.inv.halves[1,iclust,,]
    yvec = as.vector(yvec)

    ## Give the glmnet function some pre-calculated Y's and X's.
    ## fit = glmnet::glmnet(x=Xtilde[[iclust]], y=yvec,
    ##                      exclude=exclude,
    ##                      lambda=mean_lambda,
    ##                      alpha=1, intercept=FALSE, family = "gaussian")
    ## if(!is.null(coef_limit)){
    ##   res = cvxr_lasso(yvec, Xtilde[[iclust]], mean_lambda,
    ##                    exclude.from.penalty=exclude.from.penalty,
    ##                    coef_limit=coef_limit) ## Experimental feature
    ##   b=res
    ## } else {
      res = solve_lasso(x=Xtilde[[iclust]], y=yvec, lambda=mean_lambda,
                        intercept=FALSE, exclude.from.penalty=exclude.from.penalty)
    b = res$b
    ## }

    ## Obtain the coef and fitted response
    betahat = matrix(b, ncol=dimdat)
    yhat = Xa %*% betahat

    return(list(betahat=betahat, yhat=yhat))
  })

  ## Extract results
  betahats = lapply(results, function(a){a$betahat})
  yhats = lapply(results, function(a){as.matrix(a$yhat)})
  yhats = abind::abind(yhats, along=0)
  yhats = aperm(yhats, c(2,3,1)) ## This needs to by (T x dimdat x numclust)

  ## Each are lists of length |numclust|.
  return(list(beta=betahats,
              mns=yhats))
}
