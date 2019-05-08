##' Updates alpha in the covariate EM algorithm, using GLMnet as the workhorse.
## (@param beta Is a (T x p+1 x 3) array. not needed)
## (@param alpha Is a T by p+1 array. not needed)
##' @param resp Is an (T x nt x K) array.
##' @param X
##' @return The multinomial logit model coefficients. A matrix of dimension (K x
##'   (p+1)).
Mstep_alpha <- function(resp, X, numclust, lambda=0, alpha=1, iter){

  TT = nrow(X)
  p = ncol(X)

  ## Calculate the summed responsibilities
  resp.sum = t(sapply(resp, colSums)) ## (T x numclust)

  stopifnot(dim(resp)==c(TT, numclust))

  ## Fit the model, possibly with some regularization
  fit = glmnet(x=X, y=resp.sum, family="multinomial", lambda=lambda, alpha=alpha)

  ## While you're at it, calculate the fitted values (\pi) as well:
  piehat = predict(fit, newx=X, type="response")[,,1]
  stopifnot(all(dim(piehat)==c(TT,numclust)))
  ## This should be dimension (T x K)

  ## The coefficients for the multinomial logit model are (K x (p+1))
  alphahat = do.call(rbind, lapply(coef(fit), t))
  stopifnot(all(dim(alphahat)==c(numclust, (p+1))))

  return(list(pie=piehat, alpha=alphahat))
}


##' Essentially a regression
Mstep_beta <- function(resp, ylist, X){

  Xa = cbind(rep(1,nrow(X)), X)

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

  vars = lapply(1:numclust, function(iclust){

    ## Calculate all the residuals
    wt.resids = lapply(1:TT, function(tt){
      myresid = ylist[[tt]] - mn[tt,,iclust]
      sqrt(resp[[tt]][,iclust]) * myresid
    })

    ## Stack them all as rows, and take the inner product
    resid.rows = do.call(rbind, wt.resids)

    ## Calculate the covariance matrix
    resp.grandsum = sum(sapply(1:TT, function(tt){ resp[[tt]][,iclust] }))
    myvar = t(resid.rows) %*% resid.rows / resp.grandsum

    return(myvar)
  })

  ## reformat
  sigma = abind::abind(lapply(1:TT, function(tt){
    (abind::abind(vars,along=0))
  }), along=0)
  ## stopifnot(all(dim(sigma) == c(TT, numclust, dimdat, dimdat)))

  ## Return |K| of them.
  return(sigma)
}


