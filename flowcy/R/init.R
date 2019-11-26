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


##' Mean is initialized here.
init_mn <- function(ylist, numclust, TT, dimdat, warmstart=  c("none", "rough"), countslist=NULL){
  warmstart = match.arg(warmstart)
  if(warmstart == "rough"){
    return(init_mn_warmstart(ylist, numclust, countslist))
  } else if (warmstart == "none"){
    return(init_mn_naive(lapply(ylist, cbind), numclust, TT, countslist))
  } else {
    stop("warmstart option not recognized")
  }
}



##' Initialize the cluster centers (naively).
##'  @param data  A T-length list of (nt  by 3) datasets.  There should  be T of
##'   such datasets. 3 is actually \code{mulen}.
##' @param numclust Number of clusters (M).
##' @param TT total number of (training) time points.
##' @return An array of dimension (T x dimdat x M).
init_mn_naive <- function(data, numclust, TT, countslist){

  dimdat = ncol(data[[1]])

  ## Flatten the countslist
  countslist_flattened = countslist
  for(tt in 1:TT){
    thresh = quantile(countslist_flattened[[tt]], 0.6)
    countslist_flattened[[tt]][which(countslist_flattened[[tt]] >= thresh)] = thresh
  }


  ## Initialize the means by randomly sampling data from each time point.
  mulist = lapply(1:TT, function(tt){
    mydata = data[[tt]]
    nt = nrow(mydata)
    rows = sample(1:nt, numclust,
                  prob = countslist_flattened[[tt]]/sum(countslist_flattened[[tt]]))
    sampled.data = mydata[rows, , drop=FALSE]
    return(sampled.data)
  })


  ## New (T x dimdat x numclust)
  muarray = array(NA, dim=c(TT, dimdat, numclust))
  for(tt in 1:TT){
    muarray[tt,,] = mulist[[tt]]
  }

  gc()
  return(muarray)
}
##' A very rough warmstarts for covariate EM.
##' @param ylist list of data.
##' @param numclust number of clusters desired.
##' @return An array of dimension (T x dimdat x numclust).
init_mn_warmstart <- function(ylist, numclust){

  dimdat = ncol(ylist[[1]])
  TT = length(ylist)

  ## Collapse all the data
  all.y = do.call(rbind, ylist)

  ## Run k-means once on collapsed data.
  ## obj = kmeans(all.y, numclust)
  ## numclust = 5
  avg.num.rows = round(mean(sapply(ylist, nrow)))
  some.of.all.y = all.y[sample(1:nrow(all.y), avg.num.rows),]
  obj = kmeans(some.of.all.y, numclust, algorithm="MacQueen")

  ## New: warm start from truncated counts

  ## ## Plot the results (temporary)
  ## plot(some.of.all.y[,1:2], type='p',cex=0.1)
  ## points(obj$centers[,1:2], col='red', pch=16)
  ##   }

  ## Repeat it TT times and return it
  centres = array(NA, dim=c(TT, dimdat, numclust))
  for(tt in 1:TT){
    centres[tt,,] = t(obj$centers)
  }
  stopifnot(dim(centres) == c(TT, dimdat, numclust)) ## Unnecessary, but still.
  return(centres)
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


