##' Objectives.
##' @param pie matrix of mixture proportions, T by M.
##' @param mu array of dimension T by M by p.
##' @param data TT lengthed list of data
##' @param sigma array of dimension T by M by p by p.
##' @param alpha linear coefficients for regression on (log ratio of) pie.
##' @param beta linear coefficients for regression on mean.
objective_overall_cov <- function(mu, pie, sigma, data, pie_lambda=0, mean_lambda=0, alpha=0,
                                  beta=0){
  TT = length(data)
  numclust = dim(mu)[2] ## Temporary; there must be a better solution for this.
  loglikelihood_tt <- function(data, tt, mu, sigma, pie){
    dat = data[[tt]]

    ## One particle's log likelihood
    weighted.densities = sapply(1:numclust, function(iclust){

      mydat = dat
      mypie = pie[tt,iclust]
      mymu = mu[tt,iclust,]
      mysigma = as.matrix(sigma[tt,iclust,,])
      return(mypie * mvtnorm::dmvnorm(mydat,
                                      mean=mymu,
                                      sigma=mysigma,
                                      log=FALSE))

      })

    ## ## Trying to catch some errors
    ## if(any(is.nan(weighted.densities))) browser()
    ## if(rowSums(weighted.densities) < 0) browser()
    ## ## End of debug

    return(sum(log(rowSums(weighted.densities))))
  }

  ## Calculate the data likelilhood of one time point
  loglik = sapply(1:TT, function(tt){
    loglikelihood_tt(data, tt, mu, sigma, pie)
  })

  ## Return penalized
  l1norm <- function(coef){ sum(abs(coef)) }
  pen1 = pie_lambda * l1norm(as.numeric(unlist(alpha)))
  pen2 = mean_lambda * l1norm(as.numeric(unlist(beta)))
  obj = - sum(unlist(loglik)) + pen1 + pen2
  return(obj)
}
