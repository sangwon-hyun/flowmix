##' Objectives (negative log likelihood).
##'
##' @param pie matrix of mixture proportions, T by M.
##' @param mu array of dimension (T x dimdat x numclust ) ; old: T by M by p.
##' @param data TT lengthed list of data
##' @param sigma array of dimension T by M by p by p.
##' @param alpha linear coefficients for regression on (log ratio of) pie.
##' @param beta linear coefficients for regression on mean.
##' @param N Sum of the counts/biomass across all cytograms.
##'
##' @return Either the original objective value or, T-lengthed vector of
##'   objective values at each time point.
objective <- function(mu, pie, sigma,
                      ## TT, N, dimdat, numclust,
                      ylist,
                      pie_lambda=0, mean_lambda=0,
                      alpha = NULL, beta = NULL,
                      denslist_by_clust = NULL,
                      countslist = NULL,
                      each = FALSE, ## Per-cytogram likelihood.
                      sep = FALSE ## Per-particle likelihood.
                      ){

  ## Extract some things.
  TT = dim(mu)[1]
  numclust = dim(mu)[3]
  if(is.null(countslist)){
    ntlist = sapply(ylist, nrow)
    N = sum(ntlist)
  } else {
    Ntlist = sapply(countslist, sum)
    N = sum(Ntlist)
  }
  dimdat = ncol(ylist[[1]])

  ## Calculate the log likelihood
  loglik = sapply(1:TT, function(tt){
    if(is.null(denslist_by_clust)){
      return(loglikelihood_tt(ylist, tt, mu, sigma, pie, countslist, numclust))
    } else {
      return(loglikelihood_tt_precalculate(ylist, tt, denslist_by_clust, pie, countslist, numclust))
    }
  })

  ## Return penalized likelihood
  l1norm <- function(coef){ sum(abs(coef)) }

  ## Recent change: This now excludes the intercept!!!!
  pen1 = (if(!is.null(alpha)) pie_lambda * l1norm(alpha[,-1]) else 0)
  pen2 = (if(!is.null(beta)) mean_lambda * sum(sapply(beta, function(mybeta) l1norm(mybeta[-1,]))) else 0)
  obj = - 1/N * sum(unlist(loglik)) + pen1 + pen2

  ## Temporary addition: calculate per-particle, in the same format as ylist.
  if(sep){
    logliksep = sapply(1:TT, function(tt){
      return(loglikelihood_tt(ylist, tt, mu, sigma, pie, countslist, numclust,
                              sep=TRUE))
    })
  }

  ## Additional two options, temporarily used.
  if(each){
    return(unlist(loglik))
  }
  if(sep){
    return(logliksep)
  }

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
loglikelihood_tt_precalculate <- function(ylist, tt, denslist_by_clust, pie, countslist = NULL, numclust){

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
loglikelihood_tt <- function(ylist, tt, mu, sigma, pie, countslist = NULL, numclust,
                             sep=FALSE){
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
  if(sep)return(unname(log(rowSums(weighted.densities)) * counts)) ## temporary
  if(!sep) return(sum(log(rowSums(weighted.densities)) * counts))
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



##' A convenience function to calculate the CV score, or out-of-sample objective
##' function (no penalization), given new data.
##'
##' @param res flowmix object.
##' @param ylist New data.
##' @param countslist New count data.
##'
##' @return An out-of-sample, unpenalized likelihood (i.e. CV score).
objective_newdat <- function(res, ylist, countslist){

  ## ## (NOT USED NOW) The pred object.
  ## pred = predict.flowmix(res, newx = X)

  ## ## Calculate the cross-validation score.
  ## cvscore = objective(mu = pred$newmn,
  ##                     pie = pred$newpie,
  ##                     sigma = pred$sigma,
  ##                     ylist = ylist,
  ##                     countslist = countslist,
  ##                     pie_lambda = 0,
  ##                     mean_lambda = 0)

  ## Calculate the cross-validation score.
  cvscore = objective(mu = res$mn,
                      pie = res$pie,
                      sigma = res$sigma,
                      ylist = ylist,
                      countslist = countslist,
                      pie_lambda = 0,
                      mean_lambda = 0)
  return(cvscore)
}
