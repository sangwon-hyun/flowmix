##' E step for covariate EM.  Calculates the $k$'th ratio of the (pi * density)
##' of every datapoint, compared to the sum over all k=1:K. These are called
##' responsibilities (a posterieri membership probabilities).
##' @param pie matrix of component weights.
##' @param y data.
##' @param counts The number of counts per grid box. A vector of length $|I|$.
##' @import data.table
##' @return List of responsibility matrices, containing the posterior
##'   probabilities of the latent variable $Z$ (memberships to each cluster)
##'   given the parameter estimate. T-length list of (nt x dimdat)
Estep_covar <- function(mn, sigma, pie, ylist=NULL,
                        numclust,
                        denslist_by_clust=NULL,
                        first_iter = FALSE,
                        countslist = NULL## Temporary, binned approach
                        ){
  TT = length(ylist)
  ntlist = sapply(ylist, nrow)
  resp = list() ## Next up: try to /not/ do this.

  calculate_dens <- function(iclust, tt, y, mn, sigma, denslist_by_clust, first_iter,
                             bin=FALSE){## temporary

    mu <- mn[tt,,iclust] ## No problem with memory leak here.

    if(first_iter){
      dens = mvnfast::dmvn(y, mu = mu, sigma[iclust,,], log = FALSE)
    } else {
      dens = unlist(denslist_by_clust[[iclust]][[tt]])
    }
    return(dens)
  }

  ## Calculate the responsibilities at each time point, separately
  ncol.pie = ncol(pie)
  for(tt in 1:TT){
    ylist_tt = ylist[[tt]]

    ## Calculate the densities of data with respect to cluster centers
    densmat <- sapply(1:numclust,
                      calculate_dens,
                      ## Rest of arguments:
                      tt, ylist_tt, mn, sigma,
                      denslist_by_clust, first_iter,
                      bin = bin)

    ## Weight them by pie, to produce responsibilities.
    wt.densmat <- matrix(pie[tt,], nrow = ntlist[tt], ncol=ncol.pie, byrow=TRUE) * densmat
    wt.densmat = wt.densmat + 1E-10 ## Add some small number to prevent ALL zeros. Explained below..
    wt.densmat <- wt.densmat / rowSums(wt.densmat)

    ## Q: why add a small number? This is practical because it likely won't
    ## affect the end results and prevents clogging up the code. When all resp
    ## Okay, now I thnk this happens when NONE of the $n_t$ particles are likely
    ## to be in a cluster.

    resp[[tt]] <- wt.densmat
  }

  ## If |countslist| is provided, reweight the responsibilities.
  if(!is.null(countslist)){
    resp <- Map(function(myresp, mycount){ myresp * mycount },
                resp, countslist)
  }

  return(resp)
}
