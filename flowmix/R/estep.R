##' Calculates the $k$'th ratio of the (pi * density) of every datapoint,
##' compared to the sum over all clusters k=1:K. These are called responsibilities (a
##' posterieri membership probabilities).
##'
##' @param prob Matrix of component weights.
##' @param ylist Data.
##' @param mn Array of all means.
##' @param sigma (numclust x dimdat x dimdat) array.
##' @param numclust Number of clusters.
##' @param first_iter \code{TRUE} if this is the first EM iteration, which is
##'   handled separately.
##' @param denslist_by_clust Pre-calculated densities.
##' @param countslist Counts or biomass.
##'
##' @return List of responsibility matrices, containing the posterior
##'   probabilities of the latent variable $Z$ (memberships to each cluster)
##'   given the parameter estimate. T-length list of (nt x dimdat)
##'
##' @export
Estep <- function(mn, sigma, prob, ylist = NULL,
                  numclust,
                  denslist_by_clust = NULL,
                  first_iter = FALSE,
                  countslist = NULL){

  ## Setup
  TT = length(ylist)
  ntlist = sapply(ylist, nrow)
  dimdat = dim(mn)[2]

  ## Basic checks
  assertthat::assert_that(dim(mn)[1] == length(ylist))

  calculate_dens <- function(iclust, tt, y, mn, sigma, denslist_by_clust, first_iter){
    mu <- mn[tt,,iclust] ## No problem with memory leak here.
    if(first_iter){
      if(dimdat==1){
        dens = stats::dnorm(y, mu, sigma[iclust,,])
      } else {
        dens = dmvnorm_arma_fast(y, mu, sigma[iclust,,], FALSE)
      }
    } else {
      dens = unlist(denslist_by_clust[[iclust]][[tt]])
    }
    return(dens)
  }

  ## Calculate the responsibilities at each time point, separately
  ncol.prob = ncol(prob)
  resp <- lapply(1:TT, function(tt){
    ylist_tt = ylist[[tt]]

    ## Calculate the densities of data with respect to cluster centers
    densmat <- sapply(1:numclust,
                      calculate_dens,
                      ## Rest of arguments:
                      tt, ylist_tt, mn, sigma,
                      denslist_by_clust, first_iter)

    ## Weight them by prob, to produce responsibilities.
    wt.densmat <- matrix(prob[tt,], nrow = ntlist[tt], ncol = ncol.prob, byrow = TRUE) * densmat
    wt.densmat = wt.densmat + 1E-10 ## Add some small number to prevent ALL zeros.
    wt.densmat <- wt.densmat / rowSums(wt.densmat)

    ## If |countslist| is provided, reweight the responsibilities.
    if(!is.null(countslist)){
      wt.densmat = wt.densmat * countslist[[tt]]
    }
    return(wt.densmat)
  })

  return(resp)
}


