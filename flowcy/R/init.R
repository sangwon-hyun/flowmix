##' The only role is to calculate the fitted values /given/ beta.
calc_mu <- function(beta, X, dimdat, numclust){

  X.aug = cbind(1, X)

  ## We want (T x dimdat x K) fitted values, each row is a fitted mean
  fitted = array(NA, dim = c(TT, dimdat, numclust))
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
##'  @param ylist  A T-length list of (nt  by 3) datasets.  There should  be T of
##'   such datasets. 3 is actually \code{mulen}.
##' @param numclust Number of clusters (M).
##' @param TT total number of (training) time points.
##' @return An array of dimension (T x dimdat x M).
init_mn <- function(ylist, numclust, TT, dimdat, countslist = NULL){

  if(!is.null(countslist)){

    ## #### Temporary testing code ###
    ## la("flowcy")
    ## obj = make_data75()
    ## names(obj)
    ## list2env(obj, globalenv())
    ## #### End of temporary ##########

    ## NEW: initialize the means by (1) collapsing to one cytogram (2) random
    ## sampling from this distribution, after truncation,
    TT = length(ylist)
    ylist_downsampled <- lapply(1:TT, function(tt){

      y = ylist[[tt]]
      counts = countslist[[tt]]

      ## Cap the density by the median
      counts = pmax(counts, median(counts))

      ## Sample so that, in total, we get mean(nt)*30 sized sample. In the case
      ## of binned data, nt is the number of bins.
      nsize = nrow(y) / TT * 30
      y[sample(1:nrow(y), size = nsize, prob = counts/sum(counts)),, drop=FALSE]
    })

    ## Combine all and sample just |numclust| rows
    ## yy = do.call(rbind, ylist_downsampled)
    ## new_means = yy[sample(1:nrow(yy), numclust),, drop=FALSE]
    ## mulist = lapply(1:TT, function(tt){ new_means })

    ## OLD: Initialize the means by randomly sampling data from each time point.
    ## mulist = lapply(1:TT, function(tt){
    ##   y = ylist[[tt]]
    ##   nt = nrow(y)
    ##   counts = countslist[[tt]]
    ##   ## counts = countslist_flattened[[tt]]
    ##   stopifnot(length(counts) == nt)
    ##   rows = sample(1:nt, numclust,
    ##                 prob = counts / sum(counts))
    ##   sampled.data = y[rows, , drop=FALSE]
    ##   return(sampled.data)
    ## })

    ## NEW2: do a Kmeans on this.
    ## cl2 = flexclust::kcca(yy, k=numclust, family=kccaFamily("kmeans"),
    ##                       control=list(initcent="kmeanspp"))
    ## new_means = cl2@centers

    ## NEW3:
    yy = do.call(rbind, ylist_downsampled)
    new_means = yy[sample(1:nrow(yy), numclust),, drop=FALSE]

    ## Repeat TT times.
    mulist = lapply(1:TT, function(tt){ new_means })

  } else {

    mulist = lapply(1:TT, function(tt){
      y = ylist[[tt]]
      nt = nrow(y)
      rows = sample(1:nt, numclust)
      sampled.data = y[rows, , drop=FALSE]
      return(sampled.data)
    })
  }

  ## New (T x dimdat x numclust) array is created.
  muarray = array(NA, dim=c(TT, dimdat, numclust))
  for(tt in 1:TT){
    muarray[tt,,] = as.matrix(mulist[[tt]])
  }
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


