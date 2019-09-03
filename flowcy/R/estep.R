##' E step for covariate EM.  Calculates the  $k$'th ratio of the (pi * density)
##' of every  datapoint, compared to  the sum over  all k=1:K. These  are called
##' responsibilities (a posterieri membership probabilities).
##' @param pie matrix of component weights.
##' @param y data.
##' @return Responsibility matrix, containing the posterior probabilities of the
##'   latent  variable $Z$  (memberships  to each  cluster)  given the  paramter
##'   estimates. Dimension is (T X nt x K).
Estep_covar <- function(mn, sigma, pie, ylist, numclust, faster_mvn=FALSE,
                        cholspeed=FALSE,
                        eigenspeed=FALSE,
                        sigma_eig_by_dim=NULL,
                        sigma_chol_by_dim=NULL,
                        denslist_by_clust=NULL,
                        first_iter=FALSE){

  TT = length(ylist)
  ntlist = sapply(ylist, nrow)
  resp = list()
  precalculate = (eigenspeed | cholspeed)

  ## Since sigma is the same across tt, only extract it once.
  sigmalist_by_clust = lapply(1:numclust, function(iclust){
    as.matrix(sigma[1,iclust,,])
  })

  for(tt in 1:TT){
    y = ylist[[tt]]  ## These are nt rows of 3-variate measurements
    densmat = matrix(NA, nrow=ntlist[tt], ncol=numclust)
    for(iclust in 1:numclust){

      ## Gather mean at time tt
      mu = mn[tt,,iclust]

      ## Calculate the density
      if(first_iter | !precalculate){
        densmat[,iclust] = mvtnorm::dmvnorm(y,
                                            mean = mu,
                                            sigma = sigmalist_by_clust[[iclust]],
                                            log=FALSE)
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
