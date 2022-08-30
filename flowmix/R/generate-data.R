## Synopsis: Generate some complex *2d* artificial data to be driven by
## covariates. This is a direct copy of covar-artif-2d-generic.R in:
## https://github.com/robohyun66/flowcy/blob/9993d2fdc98d496a4e3295e045bd9e33cbc79ebe/main/covar-artif-2d-generic.R


##' Generic 2d or 3d synthetic data generation.
##'
##' @param p Number of covariates.
##' @param TT Number of timepoints.
##' @param fac Defaults to 1, and determines size of covariance matrix.
##' @param nt Number of particles in each time point.
##' @param dimdat Dimension of data.
##' @param seed Random seed.
##' @param prob1 Relative abundance for cluster 1 in the first half of the time
##'   points. Defaults to 3/4, in which case the first half of the time points
##'   have cluster relative abundances have a ratio of 6:1:1:0. In the later
##'   half, the fourth cluster appears, and the ratio becomes 6:1:1:1.
##'
##' @export
##'
##' @return list containing ylist, X
generate_data_generic <- function(p = 5, TT = 50, fac = 1, nt = 1000, dimdat = 2,
                                  seed = NULL,
                                  prob1 = 3/4){

  if(!is.null(seed)) set.seed(seed)
  ## dimdat = match.arg(dimdat)
  stopifnot(dimdat %in% c(2,3))

  ## Generate covariates.
  X = matrix(stats::rnorm(TT*p), ncol = p, nrow = TT)
  X[,1] = sin((1:TT)/TT * 4 * pi)
  X[,2] = c(rep(0, TT/2), rep(1, TT/2))
    ## X = scale(X)
  Xa = cbind(rep(1,TT), X)

  ## Generate coefficients. ## this should be (dimdat x p) = (2d x 2)
  beta11 = cbind(c(1,0), c(1,0)) ## Only affected by x1
  beta00 = cbind(c(0,0), c(0,0)) ## Not affected by anything
  beta10 = cbind(c(1,0), c(0,0)) ## Only the first coordinate is affected by x1
  intercept1 = c(0,0)
  intercept2 = c(1,1)
  intercept3 = c(0,5)
  intercept4 = c(5,5)
  beta1 = rbind(intercept1, beta11/3)
  beta2 = rbind(intercept2, beta10/3)
  beta3 = rbind(intercept3, beta00)
  beta4 = rbind(intercept4, beta00)
  beta1 = rbind(beta1, matrix(0, nrow = p-2, ncol = 2))
  beta2 = rbind(beta2, matrix(0, nrow = p-2, ncol = 2))
  beta3 = rbind(beta3, matrix(0, nrow = p-2, ncol = 2))
  beta4 = rbind(beta4, matrix(0, nrow = p-2, ncol = 2))

  ## Add another
  if(dimdat==3){
    beta1 = cbind(beta1, c(0,rep(0,p)))
    beta2 = cbind(beta2, c(1,rep(0,p)))
    beta3 = cbind(beta3, c(0,rep(0,p)))
    beta4 = cbind(beta4, c(5,rep(0,p)))
  }
  betalist = list(beta1, beta2, beta3, beta4)

  ## Generate the four response /components/.
  ## sigma = diag(rep(fac, 2))
  mn1 = Xa %*% beta1
  mn2 = Xa %*% beta2
  mn3 = Xa %*% beta3
  mn4 = Xa %*% beta4
  mnlist = list(mn1, mn2, mn3, mn4)

  ## Define mixture components
  pi1 = rep(prob1, TT) ## defaults to 3/4
  pi2 = rep((1-prob1)/2, TT)
  pi3 = rep((1-prob1)/2, TT)
  pi4 = c(rep(0, TT/2), rep(1/8, TT/2))
  pilist =  list(pi1, pi2, pi3, pi4)

  ## make a mixture component matrix
  pimat = do.call(cbind, pilist)
  pimat = pimat/rowSums(pimat)


  ## Define the number of points total
  ## ntbase = 1000
  ntlist = c(rep(nt, TT/2), rep(nt*9/8, TT/2))
  ntlist = apply(ntlist * cbind(pi1,pi2,pi3,pi4), 1, sum)

  ## Define the covariances
  sigma1 = diag(c(1,1,1))[1:dimdat, 1:dimdat]
  sigma2 = diag(c(10,1,1))[1:dimdat, 1:dimdat]
  sigma3 = matrix(c(3,1.5,0,
                    1.5,3,0,
                    0,0,1), ncol=3)[1:dimdat, 1:dimdat]
  sigma4 = diag(c(1,1,1))[1:dimdat, 1:dimdat]
  sigmalist = list(sigma1, sigma2, sigma3, sigma4)
  sigmalist = lapply(sigmalist, function(a) a/3*fac)
  stopifnot(all(unlist(lapply(sigmalist, dim)) == dimdat))

  ## Then, the resulting |y| is a probabistic mixture of the /components/
  datapoints = sapply(1:TT, function(tt){
    dat = get_mixture_at_timepoint(tt, ntlist[[tt]], mnlist, pilist,
                                   sigmalist = sigmalist)
  })

  ## Reformat
  ylist = lapply(datapoints, cbind)
  classlist = lapply(ylist, function(a)a[,"membership"])
  '%ni%' <- function(x, y){  return( !(x %in% y) )}
  ylist = lapply(ylist, function(a){
    a[,which(colnames(a) %ni% "membership")]
  })

  ## Return the results
  return(list(ylist = ylist,
              classlist = classlist,
              X = X,
              Xa = Xa,
              sigmalist = sigmalist,
              betalist = betalist,
              ntlist = ntlist))
}


##' Helper for generating mixture of means
##'
##' @noRd
get_mixture_at_timepoint <- function(tt, nt, mnlist, pilist, sigma=NULL, sigmalist = NULL){

  ## Check dimensions
  stopifnot(length(mnlist) == length(pilist))
  stopifnot(nrow(sigma) == ncol(mnlist[[1]]))
  dimdat = length(mnlist[[1]])

  ## Draw data from randomly chosen mean.
  prob = sapply(pilist, function(prob)prob[[tt]])
  prob = prob/sum(prob)
  stopifnot(sum(prob) == 1)
  mns = lapply(1:length(mnlist),
               function(ii){mnlist[[ii]][tt,]})
  numclust = length(prob)

  ## Samples nt memberships out of 1:numclust according to the probs in prob.
  draws = sample(1:numclust, size = nt, replace = TRUE, prob = prob)

  ## Randomly chosen means according to pi
  if(!is.null(sigmalist)){
    dat = list()
    for(iclust in 1:numclust){
      ntk = sum(draws == iclust)
      if(ntk == 0) next
      membership = rep(iclust, ntk)
      dat[[iclust]] = cbind(MASS::mvrnorm(n = ntk,
                                      mu = mns[[iclust]],
                                      Sigma = sigmalist[[iclust]]),
                        membership)
    }
    datapoints = do.call(rbind, dat)
  } else {
    ## Add noise to the means.
    means = mns[draws]
    means = do.call(rbind, means)
    datapoints = means + MASS::mvrnorm(n = nt, mu = rep(0,dimdat), Sigma = sigma)
  }

  return(datapoints)
}







##' (To add to package) add step function to the lags
##'
##' @param X Covariate matrix.
##' @param lags Number of lags.
##'
##' @return Transformed X.
##'
add_lagpar <- function(X, lags){
  ## lags = c(0,3,6,9,12)
  par = scale(X[,"par"])
  par = par - min(par)
  par = par / max(par)


  ##' Helper function to lag a vector
  lagpad <- function(x, k) {
    if (k>0) {
      return (c(rep(NA, k), x)[1 : length(x)] );
    }
    else {
      return (c(x[(-k+1) : length(x)], rep(NA, -k)));
    }
  }

  parlist = lapply(lags, function(lag) lagpad(par, lag))

  ## Make the additional columns
  dat = do.call(cbind, parlist)
  X = cbind(dat, X)
  colnames(X)[1:ncol(dat)] = c(paste0("p", 1:ncol(dat)))
  return(X)
}

##' Add a column of 0's and 1's with one changepoint, demarcated by transition
##' region crossings.
##'
##' @param X Covariates.
##' @param lat Latitude.
##' @export
add_transition <- function(X, lat){

  stopifnot(length(lat) == nrow(X))

  ## Define the crossings of the transition zone (latitude of 37).
  TT = nrow(X)
  ind = rep(0, TT)
  north = which(lat > 37)
  regions = sapply(1:3, function(ii){
    endpt = c(0,range(north),TT)[ii:(ii+1)]
    (endpt[1]+1):endpt[2]
  })

  ## Create the TF bases.
  bases0 = lapply(regions, function(reg){
    vec = rep(0, TT)
    vec[reg] = 1
    vec
  })

  ## Add them to X and return
  bases = do.call(cbind, bases0[-1])
  X = cbind(bases, X)
  colnames(X)[1:ncol(bases)] = paste0("b", 1:ncol(bases))
  return(X)
}




##' Function to create 1d simulated data.
##'
##' @param TT Number of time points.
##' @param sigma1 variance of cluster 1.
##' @param sigma2 variance of cluster 2.
##'
##' @export
generate_1d_data <- function(TT, sigma1 = NULL, sigma2 = NULL){

  ## Generate covariates.
  stopifnot(TT %% 2 == 0)
  p = 2
  X = matrix(stats::rnorm(TT*p), ncol = p, nrow = TT)
  X[,1] = sin((1:TT)/TT * 4 * pi)
  X[,2] = c(rep(0, TT/2), rep(1, TT/2))
  Xa = cbind(rep(1,TT), X)

  ## Generate coefficients. ## this should be (dimdat x p) = (2d x 2)
  offset =  3
  beta1 = c(offset + 0,1,0)
  beta2 = c(offset + 3,0,0)

  betalist = list(beta1, beta2)

  ## Generate the four response /components/.
  mn1 = Xa %*% beta1
  mn2 = Xa %*% beta2
  mnlist = list(mn1, mn2)

  ## Define mixture components
  prob1 = 3/4
  pi1 = c(rep(prob1, TT/2), rep(1E-8, TT/2))
  pi2 = c(rep((1-prob1), TT/2), rep(1-1E-9, TT/2))
  pilist =  list(pi1, pi2)

  ## make a mixture component matrix
  pimat = do.call(cbind, pilist)
  pimat = pimat/rowSums(pimat)

  ## Also define alphas
  alpha1 = c(log(3/4), 0, -8+log(3/4)) ##+10
  alpha2 = c(log(1/4), 0, -log(1/4)) ## +10
  pimat2 = cbind(exp(Xa %*% cbind(alpha1)), (exp(Xa %*% cbind(alpha2))))

  ## Note, the alpha coefficients are only identifiable up to a constant; in
  ## other words, you can plug in (alpha1 + c) and (alpha2 + c) and the
  ## resulting pi (pimat/rowSums(pimat) is the same.

  ## Define the number of points total
  nt = 10
  ntlist = c(rep(nt + (nt * 3), TT/2), rep(nt, TT/2)) %>% round()

  ## Define the covariances
  if(is.null(sigma1))sigma1 = .1
  if(is.null(sigma2))sigma2 = .1
  sigmalist = list(sigma1, sigma2)
  dimdat = 1

  ## Then, the resulting |y| is a probabistic mixture of the /components/
  datapoints = sapply(1:TT, function(tt){
    dat = get_mixture_at_timepoint(tt, ntlist[[tt]], mnlist, pilist,
                                   sigmalist = sigmalist)
  })

  ## Reformat
  ylist = lapply(datapoints, cbind)
  classlist = lapply(ylist, function(a)a[,"membership"])
  '%ni%' <- function(x, y){  return( !(x %in% y) )}
  ylist = lapply(ylist, function(a){
    a[,which(colnames(a) %ni% "membership")]
  })
  ylist = lapply(ylist, function(a) cbind(a))

  ## Return the results
  return(list(ylist = ylist,
              classlist = classlist,
              X = X,
              Xa = Xa,
              sigmalist = sigmalist,
              betalist = betalist,
              ntlist = ntlist))

}
