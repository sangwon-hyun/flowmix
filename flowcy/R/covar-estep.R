##' E step for covariate EM.  Calculates the  $k$'th ratio of the (pi * density)
##' of every  datapoint, compared to  the sum over  all k=1:K. These  are called
##' responsibilities (a posterieri membership probabilities).
##' @param pie matrix of component weights.
##' @param y data.
##' @return Responsibility matrix, containing the posterior probabilities of the
##'   latent variable $Z$  given the paramter estimates. Dimension is  (T X nt x
##'   K).
Estep_covar <- function(mn, sigma, pie, y, numclust, ntlist, iter){

  ## densities = array(NA, dim=dim(pie))
  ## wt.densities=list()
  ## wt.densmat
  ## pies = list()
  resp = list()##array(NA, dim=c(T, nt, numclust)
  for(tt in 1:TT){
    y = ylist[[tt]]  ## These are nt rows of 3-variate measurements
    densmat = matrix(NA, nrow=ntlist[tt], ncol=numclust)
    for(kk in 1:numclust){
      ## kk=2
      mu = mn[tt,,kk] ## This is a single 3-variate mean.
      sigm = sigma[[kk]]
      ## mn[tt,,]
      ## y[1]
      densmat[,kk] = mvtnorm::dmvnorm(y, mean=mu, sigma=matrix(sigm), log=FALSE)
    }
    ## if(iter==3)browser()

    piemat = matrix(pie[tt,], nrow=ntlist[tt], ncol=ncol(pie), byrow=TRUE)

    stopifnot(dim(piemat)==dim(densmat))
    wt.densmat = piemat * densmat
    wt.densmat = wt.densmat / rowSums(wt.densmat) ## make sure this is a row-wise operation
    resp[[tt]] = wt.densmat
  }
  return(resp)
}


## Design of Estep(): Calculate densities of (i, k, t) over all i=1:nt, k=1:K, t=1:T
## This is done by (1) subtracting y_{it} from $mu$
## (2) Feeding it into the 3-variate distribution $\phi()$ with the covariance $Sigma_k$.
## So, the loop should be over $t$ and $k$, and then the densities of all $i$
## should be calculated from there
