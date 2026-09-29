##' Calculates the $k$'th ratio of the (pi * density) of every datapoint,
##' compared to the sum over all clusters k=1:K. These are called responsibilities (a
##' posterieri membership probabilities).
##'
##' @param prob Matrix of component weights.
##' @param ylist Data.
##' @param mn Array of all means.
##' @param sigma (numclust x dimdat x dimdat) array.
##' @param numclust Number of clusters.
##' @param eps A small number that is added to the weighted probabilities
##'   /before/ normalizing to get the responsibilities. Defaults to 1E-20.
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
                  eps = 1E-20,
                  countslist = NULL){

  ## Setup
  TT = length(ylist)
  ntlist = sapply(ylist, nrow)
  dimdat = dim(mn)[2]

  ## Basic checks
  assertthat::assert_that(dim(mn)[1] == length(ylist))

  calculate_dens <- function(iclust, tt, y, mn, sigma, denslist_by_clust, first_iter){
   mu <- mn[tt,,iclust]
    if(first_iter){
      if(dimdat == 1){
        dens = stats::dnorm(y, mu, sd = sqrt(sigma[iclust,,])) ## make sure to use standard deviation here!
      } else {
        dens = dmvnorm_arma_fast(y, t(mu), sigma[iclust,,], FALSE)
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

    if(nrow(ylist_tt) == 0){
      return(ylist_tt)
    }

    ## Calculate the densities of data with respect to cluster centers
    densmat <- sapply(1:numclust,
                      calculate_dens,
                      ## Rest of arguments:
                      tt, ylist_tt, mn, sigma,
                      denslist_by_clust, first_iter)

    ## Weight them by prob, to produce responsibilities.
    wt.densmat <- matrix(prob[tt,], nrow = ntlist[tt], ncol = ncol.prob, byrow = TRUE) * densmat
    wt.densmat = wt.densmat + eps ## Add some small number to prevent ALL zeros.
    wt.densmat <- wt.densmat / rowSums(wt.densmat)

    ## If |countslist| is provided, reweight the responsibilities.
    if(!is.null(countslist)){
      wt.densmat = wt.densmat * countslist[[tt]]
    }
    return(wt.densmat)
  })

  return(resp)
}






Estep_new <- function(mn, sigma, prob, ylist = NULL,
                      numclust,
                      log_denslist_by_clust = NULL,
                      first_iter = FALSE,
                      eps = 1E-20,## Not used here.
                      countslist = NULL){

  ## Setup
  TT = length(ylist)
  ntlist = sapply(ylist, nrow)
  dimdat = dim(mn)[2]

  ## Basic checks
  assertthat::assert_that(dim(mn)[1] == length(ylist))

  ## We convert prob to log-prob once
  log_prob = log(prob)

  resp <- lapply(1:TT, function(tt){
    ylist_tt = ylist[[tt]]
    if(nrow(ylist_tt) == 0) return(ylist_tt)

    ## 1. Extract Log-Densities (nt x numclust)
    ## Assumes your denslist_by_clust now contains LOG densities
    log_densmat <- sapply(1:numclust, function(iclust){
      mu <- mn[tt,,iclust]
      if(first_iter){
        ## ... (Call dmvnorm_arma_fast with logd = TRUE)
        if(dimdat==1){
          ## make sure to use standard deviation here!
          log_dens = log(stats::dnorm(y, mu, sd = sqrt(sigma[iclust,,])))
        } else {
          log_dens = dmvnorm_arma_fast(ylist_tt, t(mu), sigma[iclust,,], logd = TRUE)
        }
      } else {
        unlist(log_denslist_by_clust[[iclust]][[tt]])
      }
    })

    ## 2. Calculate Log-Numerator: log(pi_k) + log(f_k)
    ## We use sweeping to add log_prob[tt,] to each row of log_densmat
    log_wt_densmat <- sweep(log_densmat, 2, log_prob[tt,], "+")

    ## 3. Log-Sum-Exp Trick for Stability
    ## Subtract the max log-value in each row to prevent exp() from blowing up or hitting zero
    row_max <- apply(log_wt_densmat, 1, max)

    # This is effectively: exp(log_wt - max) / rowSums(exp(log_wt - max))
    wt.densmat <- exp(sweep(log_wt_densmat, 1, row_max, "-"))
    wt.densmat <- wt.densmat / rowSums(wt.densmat)

    ## 4. Reweight by countslist (biomass/counts)
    if(!is.null(countslist)){
      wt.densmat = wt.densmat * countslist[[tt]]
    }

    return(wt.densmat)
  })

  return(resp)
}
