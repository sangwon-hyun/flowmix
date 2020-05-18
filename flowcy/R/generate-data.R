## Synopsis: Generate some complex *2d* artificial data to be driven by
## covariates. This is a direct copy of covar-artif-2d-generic.R in:
## https://github.com/robohyun66/flowcy/blob/9993d2fdc98d496a4e3295e045bd9e33cbc79ebe/main/covar-artif-2d-generic.R




##' Generic 2d or 3d data generation.
##'
##' @export
##'
##' @return list containing ylist, X
generate_data_generic <- function(p = 5, TT = 50, fac = 1, nt = 1000, dimdat = 2,
                                  seed = NULL){

  if(!is.null(seed)) set.seed(seed)
  ## dimdat = match.arg(dimdat)
  stopifnot(dimdat %in% c(2,3))

  ## Generate covariates.
  X = matrix(rnorm(TT*p), ncol = p, nrow = TT)
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
  pi1 = rep(3/4, TT)
  pi2 = rep(1/8, TT)
  pi3 = rep(1/8, TT)
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
get_mixture_at_timepoint <- function(tt, nt, mnlist, pilist, sigma=NULL, sigmalist = NULL){

  ## Check dimensions
  stopifnot(length(mnlist) == length(pilist))
  stopifnot(nrow(sigma) == ncol(mnlist[[1]]))
  dimdat = length(mnlist[[1]])

  ## Draw data from randomly chosen mean.
  pie = sapply(pilist, function(pie)pie[[tt]])
  pie = pie/sum(pie)
  stopifnot(sum(pie) == 1)
  mns = lapply(1:length(mnlist),
               function(ii){mnlist[[ii]][tt,]})
  numclust = length(pie)

  ## Samples nt memberships out of 1:numclust according to the probs in pie.
  draws = sample(1:numclust, size = nt, replace = TRUE, prob = pie)

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

##' Generate pseudoreal, 2-mixture 1d data from the sunlight variable.
##'
##' @param bin If \code{TRUE}, transform into binned biomass.
##'
##' @return List containing data: {ylist, X, countslist}, and true underlying
##'   coefficients and mean/probs {mnmat, pie, alpha, beta}.
generate_data_1d_pseudoreal <- function(bin = FALSE, seed=NULL, datadir="~/repos/cruisedat/export"){

  ## Setup
  TT = 100
  ntlist = c(rep(800, TT/2), rep(1000, TT/2))
  p = 3
  numclust = 2

  ## Generate covariate
  ## datadir = "~/repos/cruisedat/export"
  load(file.path(datadir, "MGL1704-hourly-only-binned.Rdata"))
  par = X[, "par"]
  par = par[!is.na(par)]
  par = ksmooth(x=1:length(par), y=par, bandwidth=5, x.points = 1:length(par))$y
  if(!is.null(seed)) set.seed(seed)
  X = cbind(scale(par[1:TT]),
            rnorm(TT),
            c(rep(0, TT/2), rep(1, TT/2)))
  colnames(X) = c("par", "noise", "cp")


  ## Beta coefficients
  beta = matrix(0, ncol = numclust, nrow = p+1)
  beta[0+1,1] = 0
  beta[1+1,1] = 3
  beta[2+1,1] = 0
  beta[0+1,2] = 10
  beta[1+1,2] = -3
  beta[2+1,2] = 0
  mnmat = cbind(1, X) %*% beta
  ## matplot(mnmat, type='l')

  ## Alpha coefficients
  alpha = matrix(0, ncol = numclust, nrow = p+1)
  alpha[0+1, 1] = 0
  alpha[1+1, 1] = 0
  alpha[2+1, 1] = 0
  alpha[3+1, 1] = 0
  alpha[0+1, 2] = -10
  alpha[1+1, 2] = 0
  alpha[2+1, 2] = 0
  alpha[3+1, 2] = 10 + log(1/4)
  colnames(alpha) = paste0("clust", 1:numclust)
  rownames(alpha) = c("intercept", "par", "noise", "cp")
  pie = exp(cbind(1,X) %*% alpha)
  pie = pie/rowSums(pie)
  ## print(round(pie,3))
  ## matplot(pie, type='l', lwd=3)

  ## Samples |nt| memberships out of (1:numclust) according to the probs in pie.
  ## Data is a probabilistic mixture from these two means, over time.
  ylist = lapply(1:TT,
                 function(tt){
                   ## print(round(pie[[tt]],3))
                   draws = sample(1:numclust,
                                  size = ntlist[tt], replace = TRUE,
                                  ## prob = c(pie[[tt]], 1-pie[[tt]]))
                                  prob = c(pie[tt,1], pie[tt,2]))
                   mns = mnmat[tt,]
                   means = mns[draws]
                   datapoints = means + rnorm(ntlist[tt], 0, 1)
                   cbind(datapoints)
                 })

  ## To bin or not!
  if(!bin){
    countslist = NULL
  } else {
    stop("Binning doesn't work yet! make_grid and bin_many_cytograms aren't written for 1d data yet.")

    ## 1. Make grid
    dat.gridsize = 30
    dat.grid = flowcy::make_grid(ylist, gridsize = dat.gridsize)

    ## 2. Bin with just counts
    obj = flowcy::bin_many_cytograms(ylist, dat.grid, mc.cores = mc.cores, verbose = TRUE)
    ybin_list = obj$ybin_list
    counts_list = obj$counts_list
    sparsecounts_list = obj$sparsecounts_list ## new

    ## 4. (NOT USED) Also obtain all the binned midpoints, in a (d^3 x 3) matrix.
    ## ybin_all = make_ybin(counts = NULL, midpoints = make_midpoints(dat.grid))
    ybin_all = obj$ybin_all ## new

    ## Assign binned data to static names |ylist| and |countslist|.
    ylist = ybin_list
    countslist = counts_list
  }

  return(list(ylist = ylist, X = X,
              countslist = countslist,
              mnmat = mnmat,
              pie = pie,
              alpha = alpha,
              beta = beta))
}


##' Create some artificial but realistic data.
##'
##' @param datadir Data directory
##' @param filename File name of the "cvres" object from the 2-76-5 experiment
##'
##' @return A list containing the generating coefficients, true means, and data
##'   (ylist, X, countslist=NULL for now).
generate_data_1d_pseudoreal_from_cv <- function(datadir, seed = NULL){

  ## Load best 1d CV result
  ## cvres = blockcv_summary(2, 76, 5, 10, nrep = 5, datadir = datadir)##, subfolder="orig")
  ## saveRDS(cvres, file=file.path("~/repos/cruisedat/export", "1d-cvres.rds"))
  ## cvres = readRDS(file=file.path(datadir, filename)) #
  ## cvres = readRDS(file=file.path("~/repos/cruisedat/export", "1d-cvres.rds"))
  cvres = readRDS(file=file.path(datadir, "1d-cvres.rds"))
  ## Save this cvres and load it from datadir
  res = cvres$bestres
  X = res$X
  numclust = res$numclust
  TT = nrow(X)
  sigmas = (cvres$pretty.sigmas)

  ## Threshold beta at 0.01
  gen_beta = do.call(cbind, res$beta)
  gen_beta %>% Matrix::Matrix(sparse=TRUE)
  small = which(abs(gen_beta[2:nrow(gen_beta),])<1E-2)
  gen_beta[2:nrow(gen_beta),][small] = 0
  gen_beta %>% Matrix::Matrix(sparse=TRUE) %>% print

  ## Threshold alpha at 0.1
  gen_alpha = t(res$alpha)
  gen_alpha %>% Matrix::Matrix(sparse=TRUE)
  small = which(abs(gen_alpha[2:nrow(gen_alpha),])<1E-1)
  gen_alpha[2:nrow(gen_alpha),][small] = 0
  gen_alpha %>% Matrix::Matrix(sparse=TRUE)


  ## Beta coefficients
  mnmat = cbind(1, X) %*% gen_beta
  matplot(mnmat)
  ## matplot(mnmat, type='l')

  ## Alpha coefficients
  pie = exp(cbind(1,X) %*% gen_alpha)
  pie = pie/rowSums(pie)
  matplot(pie, type='l', lwd=3)


  ## Load ntlist
  ntlist = rep(100, TT)

  ## Samples |nt| memberships out of (1:numclust) according to the probs in pie.
  ## Data is a probabilistic mixture from these two means, over time.
  ylist = lapply(1:TT,
                 function(tt){
                   ## print(round(pie[[tt]],3))
                   draws = sample(1:numclust,
                                  size = ntlist[tt], replace = TRUE,
                                  ## prob = c(pie[[tt]], 1-pie[[tt]]))
                                  prob = c(pie[tt,]))
                   mns = mnmat[tt,]
                   means = mns[draws]
                   noises = sapply(draws, function(iclust){ rnorm(1, 0, sigmas[iclust])})
                   datapoints = means + noises
                   cbind(datapoints)
                 })

  ## Making basic plot.
  plot_ylist(ylist, countslist=NULL, scale=FALSE)
  cols = RColorBrewer::brewer.pal(numclust, "Set2")[order(colMeans(mnmat), decreasing=TRUE)]
  matlines(mnmat, col=cols[1:4])
  for(ii in 1:numclust){
    points(mnmat[,ii], pch=16, cex=pie[,ii]*5,
           col=cols[ii])
  }

  ## Make into countslist
  ## browser()
  ## dat.gridsize = 40
  ## dat.grid = flowcy::make_grid(ylist, gridsize = dat.gridsize) ## Having this to be common among all things is important
  ## obj = flowcy::bin_many_cytograms(ylist, dat.grid, mc.cores = 8, verbose = TRUE) ## This code needs to be made into 1d data
  ## ybin_list = obj$ybin_list
  ## counts_list = obj$counts_list
  countslist = NULL

  return(list(ylist = ylist, X = X,
              countslist = countslist,
              mnmat = mnmat,
              pie = pie,
              alpha = gen_alpha,
              beta = gen_beta))
}




##' (To add to package) add step function to the lags
add_lagpar <- function(X, lags){
  ## lags = c(0,3,6,9,12)
  par = scale(X[,"par"])
  par = par - min(par)
  par = par / max(par)
  parlist = lapply(lags, function(lag)lagpad(par, lag))
  ## Make the additional columns
  dat = do.call(cbind, parlist)
  X = cbind(dat, X)
  colnames(X)[1:ncol(dat)] = c(paste0("p", 1:ncol(dat)))
  return(X)
}

##' (to add to package) Add step function bases demarcated by transition region
##' crossings.
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
