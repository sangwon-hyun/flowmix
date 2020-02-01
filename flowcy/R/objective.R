##' Objectives.
##'
##' @param pie matrix of mixture proportions, T by M.
##' @param mu array of dimension (T x dimdat x numclust ) ; old: T by M by p.
##' @param data TT lengthed list of data
##' @param sigma array of dimension T by M by p by p.
##' @param alpha linear coefficients for regression on (log ratio of) pie.
##' @param beta linear coefficients for regression on mean.
objective <- function(mu, pie, sigma,
                      ## TT, N, dimdat, numclust,
                      ylist,
                      pie_lambda=0, mean_lambda=0,
                      alpha = NULL, beta = NULL,
                      denslist_by_clust = NULL,
                      countslist = NULL){

  ## Extract some things.
  TT = dim(mu)[1]
  numclust = dim(mu)[3]
  ntlist = sapply(ylist, nrow)
  N = sum(ntlist)
  dimdat = ncol(ylist[[1]])

  ## Calculate the log likelihood
  loglik = sapply(1:TT, function(tt){
    if(is.null(denslist_by_clust)){
      return(loglikelihood_tt(ylist, tt, mu, sigma, pie, countslist))
    } else {
      return(loglikelihood_tt_precalculate(tt, denslist_by_clust, pie, countslist))
    }
  })

  ## Return penalized likelihood
  l1norm <- function(coef){ sum(abs(coef)) }

  ## Recent change: This now excludes the intercept!!!!
  pen1 = (if(!is.null(alpha)) pie_lambda * l1norm(alpha[,-1]) else 0)
  pen2 = (if(!is.null(beta)) sum(sapply(beta, function(mybeta) l1norm(mybeta[-1,]))) else 0)
  obj = - 1/N * sum(unlist(loglik)) + pen1 + pen2
  return(obj)
}

##' First helper function: Calculates one particle's log likelihood using
##' precalculated data densities.
##'
##' @param tt Time point of interest.
##' @param denslist_by_clust Pre-calculated densities.
##' @param pie Population proportions.
##'
##' @return Log likelihood.
##'
loglikelihood_tt_precalculate <- function(tt, denslist_by_clust, pie, countslist = NULL){

  ## One particle's log likelihood (weighted density)
  weighted.densities = lapply(1:numclust, function(iclust){
    return(pie[tt,iclust] * denslist_by_clust[[iclust]][[tt]])
  })
  nt = nrow(ylist[[tt]])
  counts = (if(!is.null(countslist)) countslist[[tt]] else rep(1, nt))
  return(sum(log(Reduce("+", weighted.densities)) * counts))
}

##' Second helper function: Calculates one particle's log likelihood *without*
##' any pre-calculated densities.
##'
##' @param ylist Response data.
##' @param mu Means.
##' @param pie Population proportions.
##' @param tt Time point of interest.
##' @param sigma Covariances.
##' @param countslist (Optional) Counts corresponding to \code{ylist}.
##'
##' @return Log likelihood.
##'
loglikelihood_tt <- function(ylist, tt, mu, sigma, pie, countslist = NULL){
  ### STUPID MISTAKE!!!!!!!!!!!!!! Although it's hard to see why this would affect things at iter>1

  ## One particle's log likelihood (weighted density)
  weighted.densities = sapply(1:numclust, function(iclust){
    return(pie[tt,iclust] * mvnfast::dmvn(ylist[[tt]],
                                          mu = mu[tt,,iclust],
                                          sigma = as.matrix(sigma[iclust,,]),
                                          log = FALSE))
  })
  nt = nrow(ylist[[tt]])
  counts = (if(!is.null(countslist)) countslist[[tt]] else rep(1, nt))
  return(sum(log(rowSums(weighted.densities)) * counts))
}



##' A wrapper, to get objective value at a particular time or vector of times
##' \code{tt}.
##'
##' @param tt Time points of interest (single integer value, or integer vector).
##' @param ... Rest of arguments to \code{objective()}.
##'
##' @return Objective value
objective_subset <- function(times, ...){


  ## Subset things to pass to objective()
  args = list(...)
  args$mu = (args$mu)[times,,,drop=FALSE]
  args$pie = (args$pie)[times,,drop=FALSE]
  args$ylist = (args$ylist)[times]
  args$countslist = (args$countslist)[times]
  stopifnot(all(sapply(args$ylist, nrow) == sapply(args$countslist, length)))

  ## Call the problem
  return(do.call(objective, args))
}
