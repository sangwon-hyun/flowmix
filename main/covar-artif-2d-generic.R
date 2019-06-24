## Generate some complex *2d* artificial data to be driven by covariates.

##' Helper for generating mixture of 2d means
get_2d_mixture_at_timepoint <- function(tt, nt, mnlist, pilist, sigma=NULL, sigmalist=NULL){

  ## Check dimensions
  stopifnot(length(mnlist)==length(pilist))
  stopifnot(nrow(sigma)==ncol(mnlist[[1]]))

  ## Draw data from randomly chosen mean.
  pie = sapply(pilist, function(pie)pie[[tt]])
  pie = pie/sum(pie)
  stopifnot(sum(pie)==1)
  mns = lapply(1:length(mnlist),
               function(ii){mnlist[[ii]][tt,]})
  numclust = length(pie)

  ## Samples nt memberships out of 1:numclust according to the probs in pie.
  draws = sample(1:numclust, size=nt, replace=TRUE, prob=pie)

  ## Randomly chosen means according to pi
  if(!is.null(sigmalist)){
    dat = list()
    for(iclust in 1:numclust){
      ntk = sum(draws==iclust)
      if(ntk==0) next
      dat[[iclust]] = MASS::mvrnorm(n=ntk, mu=mns[[iclust]], Sigma=sigmalist[[iclust]])
    }
   datapoints = do.call(rbind, dat)
  } else {
    ## Add noise to the means.
    means = mns[draws]
    means = do.call(rbind, means)
    datapoints = means + MASS::mvrnorm(n=nt, mu=c(0,0), Sigma=sigma)
  }

  return(datapoints)
}


TT = 50
## nt = 40
fac = 1 ## 0.3  ## 0.6   ## 1

## Here, p=2, T=50, nt = 10, dimdat = 2.

## Generate covariates.
X = matrix(rnorm(TT*p), ncol=p, nrow=TT)
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
beta1 = rbind(beta1, matrix(0, nrow=p-2, ncol=2))
beta2 = rbind(beta2, matrix(0, nrow=p-2, ncol=2))
beta3 = rbind(beta3, matrix(0, nrow=p-2, ncol=2))
beta4 = rbind(beta4, matrix(0, nrow=p-2, ncol=2))

## Generate the five response /components/.
sigma = diag(rep(fac, 2))
mn1 = Xa %*% beta1
mn2 = Xa %*% beta2
mn3 = Xa %*% beta3
mn4 = Xa %*% beta4
mnlist = list(mn1, mn2, mn3, mn4)

## Define mixture components
pi1 = rep(3/4, TT)
pi2 = rep(1/8, TT)
pi3 = rep(1/8, TT)
pi4 = c(rep(0, TT/2), rep(1/8, TT/2))
pilist = list(pi1, pi2, pi3, pi4)

## Define the number of points total
## ntbase = 1000
ntlist = c(rep(nt, TT/2), rep(nt*9/8, TT/2))
ntlist = apply(ntlist * cbind(pi1,pi2,pi3,pi4), 1, sum)

## Define the covariances
sigma1 = diag(c(1,1))
sigma2 = diag(c(10,1))
sigma3 = matrix(c(3,1.5, 1.5,3), ncol=2)
sigma4 = diag(c(1,1))
sigmalist = list(sigma1, sigma2, sigma3, sigma4)
sigmalist = lapply(sigmalist, function(a) a/3)

## Then, the resulting |y| is a probabistic mixture of the /components/
datapoints = sapply(1:TT, function(tt){
  dat = get_2d_mixture_at_timepoint(tt, ntlist[[tt]], mnlist, pilist,
                              sigmalist=sigmalist)
})

## ## Plotting to see if things are okay:
## tt = 50
## dat = get_2d_mixture_at_timepoint(tt, ntlist[[tt]], mnlist, pilist,
##                             sigmalist=sigmalist)
## plot(NA, xlim=c(-2,8), ylim=c(-2,8))
## points(dat, col='grey', pch=16)
## lapply(1:length(mnlist), function(ii){
##          mn = mnlist[[ii]]
##          points(x=mn[tt,1], y=mn[tt,2], pch=toString(ii), cex=3 )
##          points(x=mn[tt,1], y=mn[tt,2], pch=1, cex=6 )
## })



## Reformat
ylist = lapply(datapoints, cbind)
## Xrep = cbind(x1,x2)[rep(1:TT, each=nt),]

## Just defining these in case it is useful.
xlim = range(do.call(rbind,ylist)[,1])
ylim = range(do.call(rbind,ylist)[,2])
