##' Objectives.
##' @param pie matrix of mixture proportions, T by M.
##' @param mu array of dimension (T x dimdat x numclust ) ; old: T by M by p.
##' @param data TT lengthed list of data
##' @param sigma array of dimension T by M by p by p.
##' @param alpha linear coefficients for regression on (log ratio of) pie.
##' @param beta linear coefficients for regression on mean.
objective_overall_cov <- function(mu, pie, sigma,
                                  TT, dimdat, numclust,
                                  ylist,
                                  pie_lambda=0, mean_lambda=0,
                                  alpha=0, beta=0,
                                  denslist_by_clust=NULL,
                                  countslist=NULL,
                                  iter##temporary
                                  ## ylist_orig
                                  ){

  ##   ## if(iter==13) browser()
  ## ## Temporary
  ## denslist_by_clust = NULL
  ## ## End of temporary

  ## Basic checks
  stopifnot(dim(mu) == c(TT, dimdat, numclust))

  ## 1. Helper function: Calculates one particle's log likelihood using
  ## precalculated data densities.
  loglikelihood_tt_precalculate <- function(tt, denslist_by_clust, pie, countslist=NULL){

    ## One particle's log likelihood (weighted density)
    weighted.densities = lapply(1:numclust, function(iclust){
      return(pie[tt,iclust] * denslist_by_clust[[iclust]][[tt]])
    })
    if(!is.null(countslist))counts = countslist[[tt]] else counts = 1
    ## browser()
    ## log(Reduce("+", weighted.densities))
    ## weighted.densities

    ## sum(log(Reduce("+", weighted.densities)) * counts)
    ## print(sum(log(Reduce("+", weighted.densities)) * counts))

    return(sum(log(Reduce("+", weighted.densities)) * counts))
  }

  ## 2. Second helper function: Calculates one particle's log likelihood without
  ## any pre-calculated densities.
  loglikelihood_tt <- function(ylist, tt, mu, sigma, pie, countlist=NULL){
    ## One particle's log likelihood (weighted density)
    weighted.densities = sapply(1:numclust, function(iclust){

      mydat = ylist[[tt]]
      mypie = pie[tt,iclust]
      mymu = mu[tt,,iclust]
      mysigma = as.matrix(sigma[iclust,,])
      return(mypie * mvnfast::dmvn(mydat,
                                   mu=mymu,
                                   sigma=mysigma,
                                   log=FALSE))
    })
    if(!is.null(countslist)) counts = countslist[[tt]] else counts = 1

    ## ## Real data:
    ## ## One particle's log likelihood (weighted density)
    ## weighted.densities.orig = sapply(1:numclust, function(iclust){

    ##   mydat = ylist_orig[[tt]]
    ##   mypie = pie[tt,iclust]
    ##   mymu = mu[tt,,iclust]
    ##   mysigma = as.matrix(sigma[iclust,,])
    ##   return(mypie * mvnfast::dmvn(mydat,
    ##                                mu=mymu,
    ##                                sigma=mysigma,
    ##                                log=FALSE))
    ## })


    ## if(iter==8) browser()

    ## ########################
    ## ## START OF DEBUG ######
    ## ########################

    ## ## Debug: Identifying the rows that cause the problem.
    ## print("binned")
    ## print(sum(log(rowSums(weighted.densities)) * counts))
    ## print("orig")

    ## hist(log(rowSums(weighted.densities)), breaks=20)
    ## summary(log(rowSums(weighted.densities.orig)))
    ## hist((rowSums(weighted.densities.orig)), breaks=20)
    ## hist((rowSums(weighted.densities)), breaks=20)

    ## which(log(weighted.densities)>0, arr.ind=TRUE)
    ## round(weighted.densities,3)
    ## ## Okay, so all of them happen for onecluster. This must mean that clu

    ## summary(log(rowSums(weighted.densities)))

    ## round(rowSums(weighted.densities),2)
    ## which.max(rowSums(weighted.densities))
    ## weighted.densities[24,]

    ## iclust = 4
    ## numclust=2
    ## mydat = ylist[[tt]]
    ##   ## mydat = ylist_orig[[tt]]
    ## ## hist(ylist_orig[[tt]][,2])
    ## ## hist(ylist[[tt]][,2])
    ## mypie = pie[tt,iclust]
    ## mymu = mu[tt,,iclust]
    ## mysigma = as.matrix(sigma[iclust,,])
    ## mysigma = diag(c(1,1,1))
    ## a =     mvnfast::dmvn(mydat,
    ##                       mu=mymu,
    ##                       sigma=mysigma,
    ##                       log=FALSE)
    ## which.max(a)
    ## mydat = ylist[[tt]][24,]
    ## hist(a)
    ## round(mysigma,3)

    ## ## Happens often: when two things happen: the sigma of a particular cluster
    ## ## becomes very small, and when one of the coordinates of the cluster means
    ## ## becomes exactly equal to a datapoint, then the density blows up. Not sure
    ## ## /why/ this happens.

    ## summary(log(rowSums(weighted.densities.orig)))
    ## summary(log(rowSums(weighted.densities)))


    ## ## Debug: Identifying the rows that cause the problem.
    ## mu[1,,]
    ## nt = nrow(ylist[[1]])
    ## diffs = ylist[[1]] - matrix(rep(mu[1,,1],nt), nrow=nt, byrow=TRUE)
    ## summary(apply(diffs, 1, function(myrow){sqrt(sum(myrow * myrow))}))

    ## nt.orig = nrow(ylist_orig[[1]])
    ## diffs.orig = ylist_orig[[1]] - matrix(rep(mu[1,,1],nt.orig), nrow=nt.orig, byrow=TRUE)
    ## summary(apply(diffs.orig, 1, function(myrow){sqrt(sum(myrow * myrow))}))

    ## matrix(1:9, ncol=3) - c(1:3)

    ## sum(log(rowSums(weighted.densities)) * counts)
    ## print(sum(log(rowSums(weighted.densities.orig))))

    ## ######################
    ## ## END OF DEBUG ######
    ## ######################


    return(sum(log(rowSums(weighted.densities)) * counts)) ## Is this
                                                           ## multiplication
                                                           ## being done
                                                           ## correctly?

  }

  ## Calculate the log likelihood
  loglik = sapply(1:TT, function(tt){

    if(is.null(denslist_by_clust)){
      (loglikelihood_tt(ylist, tt, mu, sigma, pie, countslist))
    } else {
      ## print("binned")
      ## print(loglikelihood_tt(ylist, tt, mu, sigma, pie, countslist))
      ## print("originals")
      ## print(loglikelihood_tt(ylist_orig, tt, mu, sigma, pie))
      return(loglikelihood_tt_precalculate(tt, denslist_by_clust, pie, countslist))
    }
  })

  ## Return penalized
  l1norm <- function(coef){ sum(abs(coef)) }
  pen1 = pie_lambda * l1norm(as.numeric(unlist(alpha)))
  pen2 = mean_lambda * l1norm(as.numeric(unlist(beta)))
  obj = - sum(unlist(loglik)) + pen1 + pen2
  return(obj)
}
