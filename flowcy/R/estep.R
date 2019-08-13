##' E step for covariate EM.  Calculates the  $k$'th ratio of the (pi * density)
##' of every  datapoint, compared to  the sum over  all k=1:K. These  are called
##' responsibilities (a posterieri membership probabilities).
##' @param pie matrix of component weights.
##' @param y data.
##' @return Responsibility matrix, containing the posterior probabilities of the
##'   latent  variable $Z$  (memberships  to each  cluster)  given the  paramter
##'   estimates. Dimension is (T X nt x K).
Estep_covar <- function(mn, sigma, pie, ylist, numclust, faster_mvn=FALSE){

  TT = length(ylist)
  ntlist = sapply(ylist, nrow)
  resp = list()
  for(tt in 1:TT){
    y = ylist[[tt]]  ## These are nt rows of 3-variate measurements
    densmat = matrix(NA, nrow=ntlist[tt], ncol=numclust)
    for(kk in 1:numclust){

      ## Gather the means
      mu = mn[tt,,kk] ## This is a single 3-variate mean.
      sigm = as.matrix(sigma[tt,kk,,])
      if(faster_mvn){
        densmat[,kk] = mvnfast::dmvn(y, mu=mu, sigma=sigm, log=FALSE)##sigma=matrix(sigm)
      } else {
        densmat[,kk] = mvtnorm::dmvnorm(y, mean=mu, sigma=sigm, log=FALSE)##sigma=matrix(sigm)
      }
      if(any(is.nan(densmat[,kk]))) browser()
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
