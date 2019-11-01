##' Initialization of beta.
##' @return A (p+1 by 3 by numclust) array.
init_beta <- function(p, dimdat, numclust){
  return(array(0, dim = c(p+1, dimdat, numclust)))
}


##' Initialization of alpha.
##' @return A dimdat by p+1 array
init_alpha <- function(dimdat, p){
  return(matrix(0, ncol=(p+1), nrow=dimdat))
}


##' The only role is to calculate the fitted values /given/ beta.
calc_mu <- function(beta, X, dimdat, numclust){
  ## mu = beta %*% X ## Literally all that needs to happen, more or less.
  X.aug = cbind(rep(1, nrow(X)),X)

  ## We want (T x dimdat x K) fitted values, each row is a fitted mean
  fitted = array(NA, dim=c(TT, dimdat, numclust))
  for(iclust in 1:numclust){ ## dimdat=3
    for(idim in 1:dimdat){
      fitted[,idim,iclust] = rowSums(beta[,,idim, iclust] * X.aug) ## Elementwise operation
      ## LHS is (T x p+1)
      ## RHS is (T x p+1),
      ## So the multiple is also (T x p+1)
      ## and the rowsum is now T lengthed
    }
  }
  return(fitted)
}


calc_pie <- function(TT,numclust){

  ## For now, return a dummy of
  return(matrix(1/numclust, nrow=TT, ncol=numclust))
}


##' Initialize the cluster centers (naively).
##'  @param data  A T-length list of (nt  by 3) datasets.  There should  be T of
##'   such datasets. 3 is actually \code{mulen}.
##' @param numclust Number of clusters (M).
##' @param TT total number of (training) time points.
##' @return An array of dimension (T x dimdat x M).
init_mn <- function(data, numclust, TT){

  dimdat = ncol(data[[1]])

  ## Initialize the means by randomly sampling data from each time point.
  mulist = lapply(1:TT, function(tt){
    mydata = data[[tt]]
    nt = nrow(mydata)
    sampled.data = mydata[sample(1:nt, numclust), , drop=FALSE]
    return(sampled.data)
  })
  names(mulist) = 1:TT

  ## New (T x dimdat x numclust)
  muarray = array(NA, dim=c(TT, dimdat, numclust))
  for(tt in 1:TT){
    muarray[tt,,] = mulist[[tt]]
  }

  gc()
  return(muarray)
}



##' Initialize the covariances (naively). (TODO: TT is not needed anymore)
##' @param data The (nt by 3) datasets. There should be T of them.
##' @param numclust Number of clusters (M).
##' @param TT total number of (training) time points.
##' @return An (M by dimdat by dimdat) array containing the (dimdat by dimdat)
##'   covariances.
init_sigma <- function(data, numclust, TT, fac=1){

  ndat = nrow(data[[1]])
  pdat = ncol(data[[1]])

  sigmas = lapply(1:numclust, function(iclust){
    onesigma = diag(fac * rep(1, pdat))
    colnames(onesigma) = paste0("datcol", 1:pdat)
    rownames(onesigma) = paste0("datcol", 1:pdat)
    return(onesigma)
  })
  sigmas = abind::abind(sigmas, along=0)
  return(sigmas)
}


##' (Probably not needed) T-length list of (nt x dimdat)
init_resp <- function(ntlist, dimdat){
  lapply(ntlist, function(nt) matrix(0, nrow=nt, ncol=dimdat))
}
