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
                        bin = FALSE,
                        countslist = NULL## Temporary, binned approach
                        ){
  TT = length(ylist)
  ntlist = sapply(ylist, nrow)
  resp = list() ## Next up: try to /not/ do this.

  calculate_dens <- function(iclust, tt, y, mn, sigma, denslist_by_clust, first_iter,
                             counts=NULL, bin=FALSE){## temporary
    if(is.null(counts)) counts = rep(1, length(y))

    mu <- mn[tt,,iclust] ## No problem with memory leak here.

    if(first_iter){
      if(bin){
        dens <- counts
        inds <- which(!(counts == 0)) ## There must be a faster/better way to do this.
        dens[inds] <- mvnfast::dmvn(as.matrix(y)[inds,],
                                    mu = mu,
                                    sigma[iclust,,],
                                    log = FALSE)## * counts[inds]
      } else {
        dens = mvnfast::dmvn(y, mu = mu, sigma[iclust,,], log = FALSE)
      }
    } else {
      dens = unlist(denslist_by_clust[[iclust]][[tt]])
    }
    return(dens)
  }

  ## Calculate the responsibilities at each time point, separately
  ncol.pie = ncol(pie)
  for(tt in 1:TT){
    ylist_tt = (if(!is.null(ylist)) ylist[[tt]] else midpoints)

    ## Calculate the densities of data with respe       ct to cluster centers
    if(!bin){
      densmat <- sapply(1:numclust,
                        calculate_dens,
                        ## Rest of arguments:
                        tt, ylist_tt, mn, sigma,
                        denslist_by_clust, first_iter,
                        counts = countslist[[tt]],
                        bin = bin)
    }
    if(bin){
      denslist <- lapply(1:numclust,
                         calculate_dens,
                         ## Rest of arguments:
                         tt, ylist_tt, mn, sigma,
                         denslist_by_clust, first_iter,
                         counts = countslist[[tt]],
                         bin = bin)
      densmat = Matrix::Matrix(0, ncol=ncol.pie, nrow=ntlist[tt])
      for(iclust in 1:numclust){ ## There must be a bettr way to do it.
        densmat[,iclust] <- denslist[[iclust]]
      }
    }

    ## Weight them by pie, to produce responsibilities.
    if(bin){
      wt.densmat <- Matrix(pie[tt,], nrow = ntlist[tt], ncol=ncol.pie, byrow=TRUE) * densmat
      a = Matrix::rowSums(wt.densmat)
      ## a[which(a==0)] = 1
      wt.densmat <- wt.densmat / a
    }
    if(!bin){
      wt.densmat <- matrix(pie[tt,], nrow = ntlist[tt], ncol=ncol.pie, byrow=TRUE) * densmat
      wt.densmat <- wt.densmat / rowSums(wt.densmat)
    }
    resp[[tt]] <- wt.densmat
  }
  return(resp)
}
