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
                                  faster_mvn=FALSE,
                                  sigma_eig_by_dim=NULL,
                                  sigma_chol_by_dim=NULL,
                                  denslist_by_clust=NULL
                                  ){

  ## Basic checks
  stopifnot(dim(mu) == c(TT, dimdat, numclust))

  loglikelihood_tt <- function(data_tt, tt, mu, sigma, pie, faster_mvn, sigma_chol_by_dim){
    ## dat = data[[tt]]

    ## One particle's log likelihood
    weighted.densities = sapply(1:numclust, function(iclust){

      mydat = data_tt
      mypie = pie[tt,iclust]
      mymu = mu[tt,,iclust]
      ## mysigma = as.matrix(sigma[tt,iclust,,])
      mysigma = sigma_chol_by_dim[[iclust]]

      ## One of two ways to calculate the multivariate normal:
      if(faster_mvn){

        if(is.null(sigma_chol_by_dim)){ ## Temporary, for the first step
          isChol = FALSE
          mysigma = as.matrix(sigma[tt,iclust,,])
        } else {
          isChol = TRUE
        }

        return(mypie * mvnfast::dmvn(mydat,
                                     mu=mymu,
                                     sigma=mysigma,
                                     log=FALSE,
                                     isChol=isChol))
      } else {
        return(mypie * mvtnorm::dmvnorm(mydat,
                                        mean=mymu,
                                        sigma=mysigma,
                                        log=FALSE))
      }
      })
    return(sum(log(rowSums(weighted.densities))))
  }

  ## Temporary
  loglikelihood_tt_eigenspeed <- function(data_tt, tt, mu, pie, sigma_eig_by_dim){

    ## One particle's log likelihood
    weighted.densities = sapply(1:numclust, function(iclust){

      ## Setup
      ## mydat = data[[tt]]
      mydat = data_tt
      mypie = pie[tt,iclust]
      mymu = mu[tt,,iclust]
      mysigma_eig = sigma_eig_by_dim[[iclust]]

      ## Calculate weigthed density
      return(mypie * dmvnorm_fast(mydat,
                                  mu=mymu,
                                  sigma_eig=mysigma_eig))
    })
    return(sum(log(rowSums(weighted.densities))))
  }


  loglikelihood_tt_eigenspeed_new <- function(tt, denslist_by_clust, pie){

    ## One particle's log likelihood (weighted density)
    weighted.densities = sapply(1:numclust, function(iclust){
      return(pie[tt,iclust] * denslist_by_clust[[iclust]][[tt]])
    })
    return(sum(log(rowSums(weighted.densities))))
  }

  if(!is.null(sigma_eig_by_dim) | !is.null(sigma_chol_by_dim)){

    ## loglik = sapply(1:TT, function(tt){
    ##   loglikelihood_tt_eigenspeed(data[[tt]], tt, mu, pie, sigma_eig_by_dim)
    ## })

    loglik = sapply(1:TT, function(tt){
      loglikelihood_tt_eigenspeed_new(tt, denslist_by_clust, pie)
    })
    ## In fact, what I could do is get /all/ the densities in one swipe, without
    ## sapplying in TT. How? For a given idim, concatenate all the data (only do
    ## once), concatenate all the means (every time), concatenate the weights,
    ## calculate /all/ the densities and multiply the weights, and chop them up
    ## to store in |loglik|. This way, dmvnorm() (or its alternatives) are only
    ## called 4 times, intead of 4*TT times. This is next up!!!!
  } else {
    ## Calculate the data likelilhood of one time point
    loglik = sapply(1:TT, function(tt){
      loglikelihood_tt(data[[tt]], tt, mu, sigma, pie, faster_mvn,
                       sigma_chol_by_dim) ## Temporary addition
    })
  }
  ## End of temporary


  ## Return penalized
  l1norm <- function(coef){ sum(abs(coef)) }
  pen1 = pie_lambda * l1norm(as.numeric(unlist(alpha)))
  pen2 = mean_lambda * l1norm(as.numeric(unlist(beta)))
  obj = - sum(unlist(loglik)) + pen1 + pen2
  return(obj)
}
