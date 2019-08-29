##' E step for covariate EM.  Calculates the  $k$'th ratio of the (pi * density)
##' of every  datapoint, compared to  the sum over  all k=1:K. These  are called
##' responsibilities (a posterieri membership probabilities).
##' @param pie matrix of component weights.
##' @param y data.
##' @return Responsibility matrix, containing the posterior probabilities of the
##'   latent  variable $Z$  (memberships  to each  cluster)  given the  paramter
##'   estimates. Dimension is (T X nt x K).
Estep_covar <- function(mn, sigma, pie, ylist, numclust, faster_mvn=FALSE,
                        sigma_eig_by_dim=NULL,
                        sigma_chol_by_dim=NULL,
                        denslist_by_clust=NULL){

  TT = length(ylist)
  ntlist = sapply(ylist, nrow)
  resp = list()
  for(tt in 1:TT){
    y = ylist[[tt]]  ## These are nt rows of 3-variate measurements
    densmat = matrix(NA, nrow=ntlist[tt], ncol=numclust)
    for(iclust in 1:numclust){

      ## Gather the means
      mu = mn[tt,,iclust] ## This is a single 3-variate mean.
      ## sigm = as.matrix(sigma[tt,iclust,,]) ## This is correct
      sigm = sigma_chol_by_dim[[iclust]]
      ## By the way, I can also try not to call this many (TT) times from memory.

      if(!is.null(sigma_eig_by_dim) | !is.null(sigma_chol_by_dim)){
        mysigma_eig = sigma_eig_by_dim[[iclust]]
        ## densmat[,iclust] = dmvnorm_fast(y,
        ##                                 mu = mu,
        ##                                 sigma_eig = mysigma_eig)
        densmat[,iclust] = unlist(denslist_by_clust[[iclust]][[tt]])
      } else if(faster_mvn){
        if(is.null(sigma_chol_by_dim)){ ## Temporary, for the first step
          isChol = FALSE
          sigm = as.matrix(sigma[tt,iclust,,])
        } else {
          isChol = TRUE
        }
        densmat[,iclust] = mvnfast::dmvn(y, mu=mu, sigma=sigm, log=FALSE,
                                         isChol=isChol
                                         )##sigma=matrix(sigm)
      } else {
        densmat[,iclust] = mvtnorm::dmvnorm(y, mean=mu, sigma=sigm, log=FALSE)##sigma=matrix(sigm)
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
