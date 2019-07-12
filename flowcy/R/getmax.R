## Synopsis: Contains the functions required for getting the maximum
## regularization values.

##' Make maximum possible value of lambda regularization parameter for alpha.
##' @param  sigmalist0 |numclust| length  list of covariance matrices  (dimdat x
##'   dimdat).
##' @param alpha0 Initial intercept.
##' @param beta0 Initial intercept.
##' @param numclust Number of clusters.
##' @param ntlist Number of data points per timepoint.
##' @param X covariates (TT x dimdat)
##' @return Maximum lambda value
lambda_alpha_max <- function(alpha0, beta0, sigmalist0, numclust, X, ylist){

  ## Setup
  dimdat = nrow(beta0)
  ntlist = sapply(ylist, nrow)
  TT = length(ylist)

  ## Prep for main calcultion
  W = make_W(alpha0, numclust)
  philist = make_philist(beta0, sigmalist0, ylist, dimdat) ## TT-length list
  pielist = sapply(1:numclust, function(k){
    pi.ktilde = pifun(alpha0, k, numclust)
  })

  ## Form the matrix quantity by an entry-wise sum.
  cum.it = 1
  allmats = list()
  for(tt in 1:TT){
    for(ii in 1:nt){
      Xtrep = matrix(rep( X[tt,], numclust), ncol=dimdat, byrow=TRUE)
      numer = W %*% diag(philist[[tt]][ii,]) %*% Xtrep
      denom = sum(pielist %*% philist[[tt]][ii,])
      allmats[[cum.it]] = numer/denom
      cum.it = cum.it + 1
    }
  }

  ## Take the entrywise maximum of the sum of matrices
  sum.mat = Reduce("+", allmats)
  return(max(sum.mat))
}

##' Make maximum possible value of lambda regularization parameter for bet.
##' @param  sigmalist0 |numclust|  length list of  covariance matrices  (dimdat x
##'   dimdat).
##' @param alpha0 Initial intercept.
##' @param beta0 Initial intercept.
##' @param numclust Number of clusters.
##' @param ntlist Number of data points per timepoint.
##' @param X covariates (TT x dimdat)
##' @return Maximum lambda value
lambda_beta_max <- function(alpha0, beta0, sigmalist0, numclust, X, resplist0, ylist){

  ## Prepare things:
  ntlist = sapply(ylist, nrow)
  dimdat = nrow(beta0)
  TT = length(ylist)
  W = make_W(alpha0, numclust)
  philist = make_philist(beta0, sigmalist0, ylist, dimdat) ## TT-length list
  pielist = sapply(1:numclust, function(k){
    pi.ktilde = pifun(alpha0, k, numclust) })
  resp = resplist0

  ## Maximize this over numclust.
  sums.by.clust = list()
  for(iclust in 1:numclust){

    ## Calculate residuals from center (matrix whose rows are d-length vectors)
    resmatlist = list()
    for(tt in 1:TT){
      nt = ntlist[[tt]]
      resmatlist[[tt]] = ylist[[tt]] - matrix(rep(beta0[,iclust], nt),
                                              ncol=dimdat, byrow=TRUE)
    }

    ## Form the matrix quantity by an entry-wise sum.
    matlist2 = list()
    sigmainv = solve(sigmalist0[[iclust]])
    for(tt in 1:TT){
      matlist = lapply(1:ntlist[tt], function(ii){
        resp[[tt]][ii,iclust] * sigmainv %*% resmatlist[[tt]][ii,] %o% X[tt,]
      })
      sum.matlist = Reduce("+", matlist)
      matlist2[[tt]] = sum.matlist
    }
    sums.by.clust[[iclust]] = Reduce("+", matlist2)
  }

  ## Infin norm for matrix deriv for each cluster
  maxabs = max(sapply(sums.by.clust, function(mat)max(abs(mat))))
  return(maxabs)
}


##' Helper for getting the maximum regularization parameter.
##' @param alpha0 intercept found from standard GMM.
make_W <- function(alpha0, numclust){
  K = numclust ## replace later
  W = matrix(NA, nrow=K-1, ncol=K)
  for(ktilde in 1:(K-1)){
    pi.ktilde = pifun(alpha0, ktilde, K)
    for(k in 1:K){
      pi.k = pifun(alpha0, k, K)
      W[ktilde, k] = (if(ktilde == k){
                        pi.ktilde * (1-pi.ktilde)
                      } else {
                        W[ktilde, k] = - pi.k * pi.ktilde
                      })
    }
  }
  return(W)
}

##' Helper.
pifun <- function(inputs, k, K){
  stopifnot(length(inputs)==K-1)
  expsum =  sum(exp(inputs))
  if(k==K){
    return( 1 / (1 + expsum) )
  } else {
    return( exp(inputs[k]) / (1 + expsum))
  }
}

##' Helper for getting the maximum regularization parameter.
make_philist <- function(beta0, sigmalist, ylist, dimdat){
  K = numclust
  philist = list()
  TT = length(ylist)
  ntlist = sapply(ylist, nrow)
  for(tt in 1:TT){
    nt = ntlist[[tt]]
    phimat = matrix(NA, nrow=nt, ncol=K)
    for(kk in 1:K){
      beta0mat = matrix(rep(beta0[,kk], nt), ncol=dimdat, byrow=TRUE)
      resmat = ylist[[tt]] - beta0mat
      phimat[,kk] = mvtnorm::dmvnorm(resmat, sigma=sigmalist[[kk]]) ## nt of these
    }
    philist[[tt]] = phimat
  }
  return(philist)
}


##' Get alpha0  and beta0, which  are the corresponding  intercepts' coefficient
##' from each of the problems.
##' @param ylist List of responses
##' @param numclust Desired number of clusters.
##' @return  list containing  \code{alpha0} and  \code{beta0}, which  are fitted
##'   intercepts to the intercept-only, collapsed GMM.
get_param0 <- function(ylist, numclust){

  ## Fit model
  ylist.collapsed = do.call(rbind, ylist)
  ## if(!is.null(numpoints)){
  ##   sample.ind = sample(1:nrow(ylist.collapsed), numpoints, replace=FALSE)
  ##   ylist.collapsed = ylist.collapsed[sample.ind,]
  ## }
  res = mclust::Mclust(ylist.collapsed, numclust=numclust)

  ## Make beta0 (matrix of size d x K)
  beta0 = res$parameters$mean

  ## Make alpha0 (K-1 lengthed vector).
  resp = res$z
  pie = apply(resp, 2, sum)/sum(resp)
  alpha0 = log(pie[1:(numclust-1)]/pie[numclust])

  ## Make responsibilities, for each time.
  ntlist = sapply(ylist, nrow)
  cs = cumsum(ntlist)
  inds = Map(function(a,b){(a+1):b}, c(0,cs[-length(cs)]), cs)
  resplist0 = lapply(inds, function(ind){resp[ind,]})

  ## Make a |numclust| lengthed list of covariance matrices for each cluster
  sigmalist0 = lapply(1:numclust, function(iclust){res$parameters$variance$sigma[,,iclust]})

  return(list(alpha0=alpha0, beta0=beta0, sigmalist0=sigmalist0, resplist0=resplist0))
}
