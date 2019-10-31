##' E step for covariate EM.  Calculates the $k$'th ratio of the (pi * density)
##' of every datapoint, compared to the sum over all k=1:K. These are called
##' responsibilities (a posterieri membership probabilities).
##' @param pie matrix of component weights.
##' @param y data.
##' @return List of responsibility matrices, containing the posterior
##'   probabilities of the latent variable $Z$ (memberships to each cluster)
##'   given the parameter estimate. T-length list of (nt x dimdat)
Estep_covar <- function(mn, sigma, pie, ylist, numclust,
                        denslist_by_clust=NULL,
                        first_iter=FALSE,
                        standard_gmm=FALSE,
                        ylist_collapsed=NULL,
                        mn_collapsed=NULL
                        ){

  ## Incomplete!
  if(standard_gmm){
    assert_that(!is.null(ylist_collapsed))
    assert_that("list" %in% class(ylist_collapsed))
    TT = 1
    ntlist = nrow(ylist_collapsed[[1]])
    ylist = ylist_collapsed
    mn =  mn_collapsed
    rep = list()
    ## Also collapse means, or redefine the dimension of the means.

  } else {
    TT = length(ylist)
    ntlist = sapply(ylist, nrow)
    resp = list()
  }

  calculate_dens <- function(iclust, tt, y, mn, sigma, denslist_by_clust, first_iter){
    ## Gather mean at time tt
    mu = mn[tt,,iclust]
    ## Calculate the density
    if(first_iter){
      dens = mvnfast::dmvn(y, mu, sigma[iclust,,], log = FALSE)
    } else {
      dens = unlist(denslist_by_clust[[iclust]][[tt]])
    }
    ## if(any(is.nan(densmat[,iclust]))) browser()
    return(dens)
  }

  ## Calculate the responsibilities at each time point, separately
  for(tt in 1:TT){
    densmat <- sapply(1:numclust, calculate_dens,
                      tt, ylist[[tt]], mn, sigma, denslist_by_clust, first_iter)
    piemat <- matrix(pie[tt,], nrow=ntlist[tt], ncol=ncol(pie), byrow=TRUE)
    ## assert_that(all(dim(piemat)==dim(densmat)))
    wt.densmat <- piemat * densmat
    wt.densmat <- wt.densmat / rowSums(wt.densmat)

    ## make sure this is a row-wise operation
    resp[[tt]] = wt.densmat
  }

  return(resp)
}
