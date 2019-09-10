##' Objectives.
##' @param pie matrix of mixture proportions, T by M.
##' @param mu array of dimension (T x dimdat x numclust ) ; old: T by M by p.
##' @param data TT lengthed list of data
##' @param sigma array of dimension T by M by p by p.
##' @param alpha linear coefficients for regression on (log ratio of) pie.
##' @param beta linear coefficients for regression on mean.
objective_overall_cov <- function(mu, pie, sigma,
                                  TT, dimdat, numclust,
                                  data,
                                  pie_lambda=0, mean_lambda=0,
                                  alpha=0, beta=0,
                                  sigma_eig_by_dim=NULL,
                                  denslist_by_clust=NULL
                                  ){

  ## Basic checks
  stopifnot(dim(mu) == c(TT, dimdat, numclust))

  ## Helper function: Calculates one particle's log likelihood using
  ## precalculated data densities.
  loglikelihood_tt_precalculate <- function(tt, denslist_by_clust, pie){

    ## One particle's log likelihood (weighted density)
    weighted.densities = lapply(1:numclust, function(iclust){
      return(pie[tt,iclust] * denslist_by_clust[[iclust]][[tt]])
    })
    return(sum(log(Reduce("+", weighted.densities))))
  }

  ## Second helper function: Calculates one particle's log likelihood without
  ## any pre-calculated densities.
  loglikelihood_tt <- function(data, tt, mu, sigma, pie){

    ## One particle's log likelihood (weighted density)
    weighted.densities = sapply(1:numclust, function(iclust){

      mydat = data[[tt]]
      mypie = pie[tt,iclust]
      mymu = mu[tt,,iclust]
      mysigma = as.matrix(sigma[tt,iclust,,])
      return(mypie * mvnfast::dmvn(mydat,
                                   mu=mymu,
                                   sigma=mysigma,
                                   log=FALSE))
    })
    return(sum(log(rowSums(weighted.densities))))
  }

  loglik = sapply(1:TT, function(tt){
    if(is.null(denslist_by_clust)){
      return(loglikelihood_tt(data, tt, mu, sigma, pie))
    } else {
      return(loglikelihood_tt_precalculate(tt, denslist_by_clust, pie))
    }
  })

  ## Return penalized
  l1norm <- function(coef){ sum(abs(coef)) }
  pen1 = pie_lambda * l1norm(as.numeric(unlist(alpha)))
  pen2 = mean_lambda * l1norm(as.numeric(unlist(beta)))
  obj = - sum(unlist(loglik)) + pen1 + pen2
  return(obj)
}
