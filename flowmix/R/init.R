##' Initialize the cluster centers (naively).
##'
##' @param ylist  A T-length list of (nt  by 3) datasets.  There should  be T of
##'   such datasets. 3 is actually \code{mulen}.
##' @param numclust Number of clusters (M).
##' @param TT total number of (training) time points.
##'
##' @return An array of dimension (T x dimdat x M).
init_mn <- function(ylist, numclust, TT, dimdat, countslist = NULL){

  if(!is.null(countslist)){

    ## Initialize the means by (1) collapsing to one cytogram (2) random
    ## sampling from this distribution, after truncation,
    TT = length(ylist)
    ylist_downsampled <- lapply(1:TT, function(tt){

      y = ylist[[tt]]
      counts = countslist[[tt]]

      ## Cap the density by the median
      ## counts = pmax(counts, median(counts))

      ## Sample so that, in total, we get mean(nt)*30 sized sample. In the case
      ## of binned data, nt is the number of bins.
      nsize = pmin(nrow(y), nrow(y) / TT * 30)
      some_rows = sample(1:nrow(y), size = nsize, prob = counts/sum(counts))
      y[some_rows,, drop=FALSE]
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

    TT = length(ylist)
    ylist_downsampled <- lapply(1:TT, function(tt){
      y = ylist[[tt]]
      counts = countslist[[tt]]
      nsize = pmin(nrow(y) / TT * 30, nrow(y))
      y[sample(1:nrow(y), size = nsize),, drop=FALSE]
    })

    ## Combine all the particles
    yy = do.call(rbind, ylist_downsampled)

    ## Get K new means from these
    inds = sample(1:nrow(yy), numclust)
    new_means = yy[inds,, drop=FALSE]
    mulist = lapply(1:TT, function(tt){ new_means })

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
    ## yy = do.call(rbind, ylist_downsampled)


    ## mulist = lapply(1:TT, function(tt){
    ##   y = ylist[[tt]]
    ##   nt = nrow(y)
    ##   rows = sample(1:nt, numclust)
    ##   sampled.data = y[rows, , drop=FALSE]
    ##   return(sampled.data)
    ## })
  }

  ## New (T x dimdat x numclust) array is created.
  muarray = array(NA, dim=c(TT, dimdat, numclust))
  for(tt in 1:TT){
    muarray[tt,,] = as.matrix(mulist[[tt]])
  }
  return(muarray)
}


##' Initialize the covariances (naively).
##'
##' @param data The (nt by 3) datasets. There should be T of them.
##' @param numclust Number of clusters.
##' @param fac Value to use for the diagonal of the (dimdat x dimdat) covariance
##'   matrix.
##'
##' @return An (K x dimdat x dimdat) array containing the (dimdat by dimdat)
##'   covariances.
init_sigma <- function(data, numclust, fac = 1){

  ndat = nrow(data[[1]])
  pdat = ncol(data[[1]])
  sigmas = lapply(1:numclust, function(iclust){
    onesigma = diag(fac * rep(1, pdat))
    if(pdat==1) onesigma = as.matrix(fac)
    colnames(onesigma) = paste0("datcol", 1:pdat)
    rownames(onesigma) = paste0("datcol", 1:pdat)
    return(onesigma)
  })
  sigmas = abind::abind(sigmas, along=0)
  return(sigmas)
}
