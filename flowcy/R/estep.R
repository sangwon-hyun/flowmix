##' E step for covariate EM.  Calculates the  $k$'th ratio of the (pi * density)
##' of every  datapoint, compared to  the sum over  all k=1:K. These  are called
##' responsibilities (a posterieri membership probabilities).
##' @param pie matrix of component weights.
##' @param y data.
##' @return Responsibility matrix, containing the posterior probabilities of the
##'   latent  variable $Z$  (memberships  to each  cluster)  given the  paramter
##'   estimates. Dimension is (T X nt x K).
Estep_covar <- function(mn, sigma, pie, ylist, numclust,
                        denslist_by_clust=NULL,
                        first_iter=FALSE,
                        standard_gmm=FALSE,
                        ylist_collapsed=NULL,
                        mn_collapsed=NULL
                        ){

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


  ## Calculate the responsibilities at each time point, separately
  for(tt in 1:TT){
    y = ylist[[tt]]  ## These are nt rows of 3-variate measurements
    densmat = matrix(0, nrow=ntlist[tt], ncol=numclust)
    for(iclust in 1:numclust){

      ## Gather mean at time tt
      mu = mn[tt,,iclust]

      ## Calculate the density
      if(first_iter){
        densmat[,iclust] = mvtnorm::dmvnorm(y,
                                            mean = mu,
                                            sigma = sigma[iclust,,],
                                            log=FALSE) ## TODO: change to fastmvn::dmvn()
      } else {
        densmat[,iclust] = unlist(denslist_by_clust[[iclust]][[tt]])
      }
      if(any(is.nan(densmat[,iclust]))) browser()
    }
    piemat = matrix(pie[tt,], nrow=ntlist[tt], ncol=ncol(pie), byrow=TRUE)

    stopifnot(dim(piemat)==dim(densmat))

    wt.densmat = piemat * densmat
    wt.densmat = wt.densmat / rowSums(wt.densmat)

    ## make sure this is a row-wise operation
    resp[[tt]] = wt.densmat
  }
  return(resp)
}
