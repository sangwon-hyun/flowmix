##' Generate pseudoreal, 2-mixture 1d data from the sunlight variable.
##'
##' @param bin If \code{TRUE}, transform into binned biomass.
##' @param nt Number of points in the first cluster.
##' @param df Degrees of freedom in t distribution.
##' @param noisetype Noise type. Default is Gaussian.
##' @param dat.gridsize Grid size for binning.
##' @param datadir Directory that contains the data.
##' @param beta_par Parameter.
##' @param gap Size of gap between the two clusters. Defaults to 3.
##' @param skew_alpha skewness of skew-normal distribution that replaces
##'   N(0,1). This distribution is further shifted and scaled to have mean 0 and
##'   variance 1.
##' @param df degrees of freedom for t-distribution that replaces N(0,1). This
##'   distribution is further scaled to have mean 0 and variance 1.
##' @param no_changepoint If TRUE, remove the changepoint covariate
##'   entirely. Defaults to FALSE.
##' @param gradual_change If TRUE, the cluster probabilities transition to
##'   low/high from high/low (and vice versa) in the middle. Defaults to FALSE.
##'
##' @return List containing data: {ylist, X, countslist}, and true underlying
##'   coefficients and mean/probs {mnmat, prob, alpha, beta}.
##'
##' @export
generate_data_1d_pseudoreal <- function(bin = FALSE, seed=NULL, datadir="~/repos/cruisedat/export",
                                        nt = 1000,
                                        beta_par = 0.5,
                                        p = 3, ## Number of total covariates
                                        dat.gridsize = 30,
                                        noisetype = c("gaussian", "heavytail", "skewed"),
                                        df = NULL,
                                        skew_alpha = NULL,
                                        gap = 3,
                                        TT = 100,
                                        no_changepoint = FALSE,
                                        gradual_change = FALSE){

  ## Setup and basic checks
  assertthat::assert_that(nt %% 5 ==0)
  ntlist = c(rep(0.8 * nt, TT/2), rep(nt, TT/2))
  numclust = 2
  stopifnot(p >= 3)
  noisetype = match.arg(noisetype)
  if (noisetype == "heavytail"){ assertthat::assert_that(!is.null(df))  }
  if (noisetype == "skewed"){ assertthat::assert_that(!is.null(skew_alpha))  }

  ## Generate covariate
  ## datadir = "~/repos/cruisedat/export"
  load(file.path(datadir, "MGL1704-hourly-only-binned.Rdata"))
  par = X[, "par"]
  par = par[!is.na(par)]
  par = stats::ksmooth(x = 1:length(par), y = par, bandwidth = 5, x.points = 1:length(par))$y
  ## sst = X[,"sst"] ## In the future, maybe add another covariate

  if(!is.null(seed)) set.seed(seed)
  Xrest = do.call(cbind, lapply(1:(p-2), function(ii) stats::rnorm(TT)) )
  X = cbind(scale(par[1:TT]),
            c(rep(0, TT/2), rep(1, TT/2)),
            Xrest)## p-2 columns
  colnames(X) = c("par", "cp", paste0("noise", 1:(p-2)))

  ## Beta coefficients
  beta = matrix(0, ncol = numclust, nrow = p+1)
  beta[0+1,1] = 0
  beta[1+1,1] = beta_par
  beta[0+1,2] = gap
  beta[1+1,2] = -beta_par
  colnames(beta) = paste0("clust", 1:numclust)
  rownames(beta) = c("intercept", "par", "cp", paste0("noise", 1:(p-2)))

  ## alpha coefficients
  alpha = matrix(0, ncol = numclust, nrow = p+1)
  alpha[0+1, 2] = -10
  alpha[2+1, 2] = 10 + log(1/4)

  colnames(alpha) = paste0("clust", 1:numclust)
  rownames(alpha) = c("intercept", "par", "cp", paste0("noise", 1:(p-2)))

  ## Generate means and probabilities
  if(no_changepoint){
    ## Only designed for TT=100
    beta = beta[-3,]
    alpha = alpha[-3,]
    alpha[1, 2] = -1.4
    X = X[,-2]
  }

  if(gradual_change){
    X[,2] = c(1:(TT/2), (TT/2):1)
    alpha[0+1, 2] = -1.4
    ## a = 0.05
    a = 0.05 * 100 / TT
    alpha[2+1, 2] = a
  }

  mnmat = cbind(1, X) %*% beta
  prob = exp(cbind(1,X) %*% alpha)
  prob = prob/rowSums(prob)

  ## Samples |nt| memberships out of (1:numclust) according to the probs in prob.
  ## Data is a probabilistic mixture from these two means, over time.
  ylist = lapply(1:TT,
                 function(tt){
                   draws = sample(1:numclust,
                                  size = ntlist[tt], replace = TRUE,
                                  ## prob = c(prob[[tt]], 1-prob[[tt]]))
                                  prob = c(prob[tt,1], prob[tt,2]))
                   mns = mnmat[tt,]
                   means = mns[draws]

                   ## Add noise to obtain data points.
                   if(noisetype == "gaussian"){
                     noise = stats::rnorm(ntlist[tt], 0, 1)
                   } else if (noisetype == "heavytail"){
                     assertthat::assert_that(df >= 3)
                     variance = df / (df - 2)
                     noise = stats::rt(ntlist[tt], df = df) / sqrt(variance)
                   } else if (noisetype == "skewed"){
                     omega = sqrt(1/(1 - 2 * (1/pi) * skew_alpha^2 / (1 + skew_alpha^2)))
                     mn_shift = omega * skew_alpha * (1 / sqrt(1+skew_alpha^2)) * sqrt(2/pi)
                     noise = sn::rsn(ntlist[tt], xi = 0, omega = omega, alpha = skew_alpha) - mn_shift
                   } else {
                     stop("not written yet")
                   }
                   datapoints = means + noise
                   cbind(datapoints)
                 })

  ## To bin or not!
  if(!bin){
    countslist = NULL
  } else {
    ## stop("Binning doesn't work yet! make_grid and bin_many_cytograms aren't written for 1d data yet.")

    ## 1. Make grid
    dat.grid = flowmix::make_grid(ylist, gridsize = dat.gridsize)

    ## 2. Bin with just counts
    mc.cores = 4
    obj = flowmix::bin_many_cytograms(ylist, dat.grid, mc.cores = mc.cores, verbose = TRUE)
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


  ## Other things about the true generating model
  sigma = array(NA, dim=c(2,1,1))
  sigma[1,1,1] = sigma[2,1,1] = 1
  mn = array(NA,dim=c(100,1,2))
  mn[,1,] = mnmat
  numclust=2

  return(list(ylist = ylist, X = X,
              countslist = countslist,
              ## The true generating model:
              mnmat = mnmat,
              prob = prob,
              mn = mn,
              numclust = numclust,
              alpha = alpha,
              beta = beta,
              sigma = sigma))
}
