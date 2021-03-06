##' Main function for our method. Repeats the EM algorithm with |nrep| restarts
##' (5 by default).
##'
##' @param ... Arguments for \code{flowmix_once()}.
##' @param nrep Number of restarts.
##'
##' @return The |flowmix| class object that had the best likelihood.
##' @export
flowmix <- function(..., nrep = 5){

  ## Don't do many restarts if warmstart-able mean exists
  dots <- list(...)
  if(!is.null(dots$mn)) nrep = 1

  ## Do |nrep| restarts.
  reslist = list()
  for(itrial in 1:nrep){
    reslist[[itrial]] = flowmix_once(...)
  }

  ## Pick the best one and return
  objlist = lapply(reslist, function(res){ res$obj[-1]})
  ii = which.min(sapply(objlist, min))
  return(reslist[[ii]])
}

##' Main function for running the EM algorithm once.
##'
##' @param ylist T-length list each containing response matrices of size (nt x
##'   3), which contains coordinates of the 3-variate particles, organized over
##'   time (T) and with (nt) particles at every time.
##' @param countslist Multiplicity for particles in \code{ylist}.
##' @param numclust Number of clusters
##' @param X Matrix of size (T x p+1)
##' @param prob.list (T by K)
##' @param mean_lambda lambda for lasso for the mean.
##' @param prob_lambda lambda for lasso for probabilities.
##' @param tol_em Relative tolerance for EM convergence. Defaults to 1E-4.
##' @param zero_stabilize Defaults to FALSE. If TRUE, the EM is only run until
##'   the zero pattern in the coefficients stabilize.
##'
##' @return List containing fitted parameters and means and mixture weights,
##'   across algorithm iterations. \code{beta} is a list of (p+1 x dimdat)
##'   arrays. \code{alpha} is a (numclust x (p+1)) array.
##'
##' @export
flowmix_once <- function(ylist, X,
                         countslist = NULL,
                         numclust, niter = 1000,
                         mn = NULL, prob_lambda,
                         mean_lambda, verbose = FALSE,
                         sigma_fac = 1, tol_em = 1E-4,
                         maxdev = NULL,
                         countslist_overwrite = NULL,
                         zero_stabilize  = FALSE,
                         init_mn_flatten = FALSE,
                         ## beta Mstep (CVXR) settings
                         mstep_cvxr_ecos_thresh = 1E-8,
                         mstep_cvxr_scs_eps = 1E-5,
                         zerothresh = 1E-6,
                         ## beta Mstep (ADMM) settings
                         admm = TRUE,
                         admm_rho = 0.01,
                         admm_err_rel = 1E-3,
                         admm_err_abs = 1E-4,
                         ## beta M step (Locally Adaptive ADMM) settings
                         admm_local_adapt = TRUE,
                         admm_local_adapt_niter = 10,
                         admm_niter = (if(admm_local_adapt)1E3 else 1E4),
                         CVXR =FALSE ## temporary
                         ){

  ## Basic checks
  if(!is.null(maxdev)){
    assertthat::assert_that(maxdev!=0)
  } else {
    maxdev = 1E10 ## Some large number
  }
  ## assert_that(!(is.data.frame(ylist[[1]])))
  assertthat::assert_that(!(is.data.frame(X)))
  assertthat::assert_that(sum(is.na(X)) == 0)
  assertthat::assert_that(length(ylist) == nrow(X))
  ## assertthat::assert_that(prob_lambda > 0)
  assertthat::assert_that(numclust > 1)
  assertthat::assert_that(niter > 1)

  ## assertthat::assert_that(all(sapply(ylist, nrow) == sapply(countslist, length)))

  ## Setup
  TT = length(ylist)
  dimdat = ncol(ylist[[1]])
  p = ncol(X)
  if(is.null(mn)) mn = init_mn(ylist, numclust, TT, dimdat, countslist)
  ntlist = sapply(ylist, nrow)
  N = sum(ntlist)

  ## Initialize some objects
  prob = matrix(1/numclust, nrow = TT, ncol = numclust) ## Initialize to all 1/K.
  denslist_by_clust <- NULL
  objectives = c(+1E20, rep(NA, niter-1))
  sigma = init_sigma(ylist, numclust, sigma_fac) ## (T x numclust x dimdat x dimdat)
  sigma_eig_by_clust = NULL
  zero.betas = zero.alphas = list()
  admm_niters = list()

  ## Warm startable variables
  betas = NULL
  Zs = NULL
  wvecs = NULL
  uws = NULL
  Uzs = NULL

  ## New ADMM parameters.
  Zs = NULL
  Ws = NULL
  Us = NULL

  ## The least elegant solution I can think of.. used only for blocked cv
  if(!is.null(countslist_overwrite)) countslist = countslist_overwrite
  if(!is.null(countslist)) check_trim(ylist, countslist)

  start.time = Sys.time()
  for(iter in 2:niter){
    if(verbose){
      print_progress(iter-1, niter-1, "EM iterations.", start.time = start.time)
    }
    resp <- Estep(mn, sigma, prob, ylist = ylist, numclust = numclust,
                  denslist_by_clust = denslist_by_clust,
                  first_iter = (iter == 2), countslist = countslist)

    ## M step (three parts)
    ## 1. Alpha
    res.alpha = Mstep_alpha(resp, X, numclust, lambda = prob_lambda,
                            zerothresh = zerothresh)
    prob = res.alpha$prob
    alpha = res.alpha$alpha
    rm(res.alpha)

    ## 2. Beta

    ## temporary
    if(CVXR){
    res.beta = Mstep_beta(resp, ylist, X,
                          mean_lambda = mean_lambda,
                          first_iter = (iter == 2),
                          sigma_eig_by_clust = sigma_eig_by_clust,
                          sigma = sigma, maxdev = maxdev)
    } else {
    res.beta = Mstep_beta_admm(resp, ylist, X,
                               mean_lambda = mean_lambda,
                               first_iter = (iter == 2),
                               sigma_eig_by_clust = sigma_eig_by_clust,
                               sigma = sigma, maxdev = maxdev, rho = admm_rho,
                               betas = betas,
                               Zs = Zs,
                               Ws = Ws,
                               Us = Us,
                               err_rel = admm_err_rel,
                               err_abs = admm_err_abs,
                               niter = admm_niter,
                               local_adapt = admm_local_adapt,
                               local_adapt_niter = admm_local_adapt_niter)
  }

    admm_niters[[iter]] = unlist(res.beta$admm_niters)

    ## Harvest means
    mn = res.beta$mns
    betas = beta = res.beta$beta

    ## Harvest other things for next iteration's ADMM.
    Zs = res.beta$Zs
    Ws = res.beta$Ws
    Us = res.beta$Us
    ## rm(res.beta)

    ## Check if the number of zeros in the alphas and betas have stabilized.
    zero.betas[[iter]] = lapply(beta, function(mybeta) which(mybeta==0))
    zero.alphas[[iter]] = which(alpha == 0)
    if(zero_stabilize & iter >= 30){ ## If 5 is to low, try 10 instead of 5.
      if(check_zero_stabilize(zero.betas, zero.alphas, iter)) break
    }

    ## 3. Sigma
    sigma = Mstep_sigma(resp, ylist, mn, numclust)

    ## 3. (Continue) Decompose the sigmas.
    sigma_eig_by_clust <- eigendecomp_sigma_array(sigma)
    denslist_by_clust <- make_denslist_eigen(ylist, mn, TT, dimdat, numclust,
                                             sigma_eig_by_clust,
                                             countslist)

    ## Calculate the objectives
    objectives[iter] = objective(mn, prob, sigma, ylist,
                                 prob_lambda = prob_lambda,
                                 mean_lambda = mean_lambda,
                                 alpha = alpha, beta = beta,
                                 denslist_by_clust = denslist_by_clust,
                                 countslist = countslist)

    ## Check convergence
    ## if(iter > 10){ ## don't stop super early. ## We might not need this.
      if(check_converge_rel(objectives[iter-1],
                            objectives[iter],
                            tol = tol_em)) break
    ## }
    ## if(objectives[iter] > objectives[iter-1] * 1.01 ) break # Additional stopping
                                        ## of the likelihood
                                        ## increasing more
                                        ## than 1%.
  }

  ## Measure time
  lapsetime = difftime(Sys.time(), start.time, units = "secs")
  time_per_iter = lapsetime / (iter-1)

  ## Also calculate per-cytogram likelihoods (NOT divided by nt)
  loglikelihoods = objective(mn, prob, sigma, ylist,
                             prob_lambda = prob_lambda,
                             mean_lambda = mean_lambda,
                             alpha = alpha, beta = beta,
                             denslist_by_clust = denslist_by_clust,
                             countslist = countslist,
                             each = TRUE)
  ## loglikelihoods_particle = objective(mn, prob, sigma, ylist,
  ##                            prob_lambda = prob_lambda,
  ##                            mean_lambda = mean_lambda,
  ##                            alpha = alpha, beta = beta,
  ##                            denslist_by_clust = denslist_by_clust,
  ##                            countslist = countslist,
  ##                            each = FALSE,
  ##                            sep=TRUE)

  ## Also reformat the coefficients
  obj <- reformat_coef(alpha, beta, numclust, dimdat, X)
  alpha = obj$alpha
  beta = obj$beta

  return(structure(list(alpha = alpha,
                        beta = beta,
                        mn = mn,
                        prob = prob,
                        sigma = sigma,
                        ## denslist_by_clust = denslist_by_clust,
                        objectives = objectives[2:iter],
                        final.iter = iter,
                        time_per_iter = time_per_iter,
                        total_time = lapsetime,
                        loglikelihoods = loglikelihoods,
                        ## loglikelihoods_particle = loglikelihoods_particle,
                        ## Above is output, below are data/algorithm settings.
                        dimdat = dimdat,
                        TT = TT,
                        N = N,
                        p = p,
                        numclust = numclust,
                        X = X,
                        prob_lambda = prob_lambda,
                        mean_lambda = mean_lambda,
                        maxdev=maxdev,
                        niter = niter,
                        admm_niters = admm_niters
                        ), class = "flowmix"))
}

## ## Some tests to add
## ## object is the result of having run flowmix() or flowmix_once().
## check_size <- function(obj){
##   assert_that(check_beta_size(res$beta, p, dimdat, numclust))
##   assert_that(check_alpha_size(res$alpha, p, dimdat))
## }

## check_beta_size <- function(beta, p, dimdat, numclust){
##   all.equal(dim(beta), c(p+1, dimdat, numclust))
## }
## check_alpha_size <- function(alpha, p, dimdat){
##   all.equal(dim(alpha), c(dimdat, p+1))
## }



##' Prediction: Given new X's, generate a set of means and probs (and return the
##' same Sigma).
##'
##' @param res Object returned from covariate EM flowmix().
##' @param newx New covariate.
##'
##' @return List containing mean, prob, and sigma.
##'
##' @export
##'
predict.flowmix <- function(res, newx = NULL){

  ## ## Check the dimensions
  ## stopifnot(ncol(new.x) == ncol(res$X))
  ## newx = X[1,,drop=FALSE]
  if(is.null(newx)){
    newx = res$X
  }

  ## Check if the variable names are the same.
  cnames = res$X %>% colnames()
  cnames_new = newx %>% colnames()
  stopifnot(all(cnames == cnames_new))

  ## Augment it with a dummy variable 1
  if(nrow(newx)>1){
    newx.a = cbind(rep(1, nrow(newx)), newx)
  } else {
    newx.a = c(1, newx)
  }

  TT = nrow(newx) ## This used to be nrow(X)..
  numclust = res$numclust
  dimdat = res$dimdat
  if(is.null(dimdat)) dimdat = res %>%.$mn %>% dim() %>% .[2] ## for back=compatibility

  ## Predict the means (manually).
  newmn = lapply(1:numclust, function(iclust){
    newx.a %*% res$beta[[iclust]]
  })
  newmn_array = array(NA, dim=c(TT, dimdat, numclust))
  for(iclust in 1:numclust){ newmn_array[,,iclust] = newmn[[iclust]] }

  ## Predict the probs.
  ## newprob = predict(res$alpha.fit, newx=newx, type='response')[,,1]

  probhatmat = as.matrix(exp(cbind(1,newx) %*% t(res$alpha)))
  newprob = probhatmat / rowSums(probhatmat)
  ## predict(fit, newx=X, type="response")[,,1]
  stopifnot(all(dim(newprob) == c(TT,numclust)))
  stopifnot(all(newprob >= 0))

  ## Return all three things
  return(list(mn = newmn_array,
              prob = newprob,
              pie = newprob, ## Just a copy of prob, for back-compatibility
              alpha = res$alpha,
              beta = res$beta,
              sigma = res$sigma,
              TT = res$TT,
              N = res$N,
              numclust = res$numclust,
              X = newx))
}



##' Helper for making list of densities. Returns list by cluster then time
##' e.g. access by \code{denslist_by_clust[[iclust]][[tt]]}
##'
##' @param ylist T-length list each containing response matrices of size (nt x
##'   3), which contains coordinates of the 3-variate particles, organized over
##'   time (T) and with (nt) particles at every time.
##' @param mu (T x dimdat x numclust) array.
##' @param dimdat dimension of data.
##' @param numclust number of clusters.
##' @param TT number of time points
##' @param sigma_eig_by_clust Result of running
##'   \code{eigendecomp_sigma_array(sigma.list[[iter]])}.
##'
##' @return numclust-lengthed list of TT-lengthed.
##'
make_denslist_eigen <- function(ylist, mu,
                                TT, dimdat, numclust,
                                sigma_eig_by_clust,
                                countslist){ ## Temporary

  ## Basic checks
  assertthat::assert_that(!is.null(sigma_eig_by_clust))

  ## Calculate densities (note to self: nested for loop poses no problems)
  lapply(1:numclust, function(iclust){
    mysigma_eig <- sigma_eig_by_clust[[iclust]]
      lapply(1:TT, function(tt){
        ## return(dmvnorm_fast(ylist[[tt]],
        ##                     mu[tt,,iclust],
        ##                     sigma_eig=mysigma_eig))
        mn = mu[tt,,iclust]
        sgm = mysigma_eig$sigma
        if(dimdat == 1){
          mn = as.matrix(mn)
          sgm = sgm %>% as.matrix()
        }
        return(dmvnorm_arma_fast(ylist[[tt]],
                                 mn,
                                 sgm))
    })
  })
}





##' Functions to check convergence.
check_converge_rel <- function(old, new, tol=1E-6){ return(abs((old-new)/old) < tol )  }
