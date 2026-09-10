##' Main function for our method. Repeats the EM algorithm with |nrep| restarts
##' (5 by default).
##'
##' @param ... Arguments for \code{flowmix_once()}.
##' @param nrep Number of restarts.
##'
##' @return The |flowmix| class object that had the best likelihood.
##' @export
flowmix <- function(..., nrep = 5){

  dots <- list(...)
  if("verbose" %in% names(dots)){
    if(dots$verbose){
      cat("EM will restart", nrep, "times", fill=TRUE)
    }
  }

  ## Don't do many restarts if warmstart-able mean is provided
  if(!is.null(dots$mn)) nrep = 1
  if(!is.null(dots$seed)) stop("Can't provide seed for flowmix()! Only for flowmix_once().")

  ## Do |nrep| restarts.
  reslist = list()
  for(irep in 1:nrep){
    if("verbose" %in% names(dots)){
      if(dots$verbose){ cat("EM restart:", irep, fill=TRUE) }
    }
    reslist[[irep]] = flowmix_once(...)
    if("verbose" %in% names(dots)){
      if(dots$verbose){
        cat(fill=TRUE)
      }
    }
  }

  ## Pick the best one and return
  objlist = lapply(reslist, function(res){ res$obj[-1]})
  ii = which.min(sapply(objlist, min))
  final_model = reslist[[ii]]

  ## Also save /all/ the objectives
  final_model$all_objectives =
    lapply(1:nrep, function(irep){
        one_model = reslist[[irep]]
        data.frame(
          irep = irep, 
          iter = seq_along(one_model$objectives), 
          objective = one_model$objectives
        )
    }) %>% dplyr::bind_rows()

  return(final_model)
}

##' Main function for running the EM algorithm once.
##'
##' @param ylist T-length list each containing response matrices of size (nt x
##'   d), which contains coordinates of the d-variate particles, organized over
##'   time (T) and with (nt) particles at every time point.
##' @param countslist Multiplicity for particles in \code{ylist}.
##' @param numclust Number of clusters
##' @param X Matrix of size (T x p+1)
##' @param mean_lambda lambda for lasso for the mean.
##' @param prob_lambda lambda for lasso for probabilities.
##' @param n_princomps Number of principal components to use if principal components 
##' are to be used as regressors, rather than the original covariates (columns of X). 
##' If \code{NULL}, the original covariates will be used as regressors.
##' @param center Whether to center the columns of \code{X} before applying principal components 
##' analysis (PCA).
##' @param scale Whether to scale the columns of \code{X} before applying PCA.
##' @param n_hidden_nodes Number of hidden nodes to use in the random weight neural network (NN). If 
##' \code{NULL}, no NN is applied as the original covariates are used as regressors.
##' @param hidden_node_dist Distribution from which the hidden layer weights are randomly drawn. 
##' Currently supported options are Uniform(-\code{scale}, \code{scale}) (\code{"unif"}) and 
##' Normal(0, \code{scale}).
##' @param hidden_node_scale Scale parameter of \code{hidden_node_dist}.
##' @param hidden_node_seed RNG seed for drawing hidden layer weights.
##' @param tol_em Relative tolerance for EM convergence. Defaults to 1E-4.
##' @param zero_stabilize Defaults to FALSE. If TRUE, the EM is only run until
##'   the pattern of zeros in the coefficients stabilizes over EM iterations.
##' @param seed Seven integers that is called directly before \code{init_mn()}
##'   so it can be assigned to \code{.Random.seed} for setting the random state
##'   for the initial mean generation.
##' @param niter Number of EM iterations.
##' @param mn Initial means to use. Defaults to NULL.
##' @param verbose TRUE for loudness (e.g. printing EM iterations).
##' @param maxdev Radius for maximum deviation of cluster means over time.
##' @param admm_rho Step size for ADMM
##' @param admm_err_rel Relative error threshold for stopping ADMM.
##' @param admm_err_abs Absolute error threshold for stopping ADMM.
##' @param admm_local_adapt_niter Absolute error threshold for stopping ADMM.
##' @param admm_niter Number of ADMM iterations.
##' @param admm_local_adapt TRUE if locally adaptive ADMM (LA-ADMM) is to be used. If
##'   so, \code{admm_niter} becomes the inner number of iterations, and
##'   \code{admm_local_adapt_niter} becomes the number of outer iterations.
##' @param admm_local_adapt_niter Number of inner iterations in LA ADMM.
##' @param flatX_thresh Threshold for detecting if any covariates are flat (low
##'   variance). These flat coefficients will have be set to zero and excluded
##'   from estimation altogether.
##' @param sigma_fac Defaults to 1, and governs how big the initial covariance
##'   matrices should be in size.
##' @param countslist_overwrite Mainly used by \code{cv.flowmix()}; not to be
##'   modified by user.
##' @param zerothresh Alpha coefficient values below \code{zerothresh} are set
##'   to zero.
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
                         mean_lambda,
                         n_princomps = NULL, 
                         center = TRUE, 
                         scale = TRUE,
                         n_hidden_nodes = NULL,
                         hidden_node_dist = c("unif", "norm"), 
                         hidden_node_scale = 0.5, 
                         hidden_node_seed = NULL,
                         verbose = FALSE,
                         sigma_fac = 1, tol_em = 1E-4,
                         maxdev = NULL,
                         countslist_overwrite = NULL,
                         zero_stabilize  = FALSE,
                         zerothresh = 1E-6,
                         ## beta Mstep (ADMM) settings
                         admm_rho = 0.01,
                         admm_err_rel = 1E-3,
                         admm_err_abs = 1E-4,
                         ## beta M step (Locally Adaptive ADMM) settings
                         admm_local_adapt = TRUE,
                         admm_local_adapt_niter = 10,
                         admm_niter = (if(admm_local_adapt)1E3 else 1E4),
                         seed = NULL,
                         flatX_thresh = 1e-5
                         ){


  . = NULL ## Fixing check()

  ## Capture all arguments once
  call <- sys.call();
  call[[1]] <- as.name('list');
  args <- eval.parent(call)

  ## Basic checks
  if(!is.null(maxdev)){
    assertthat::assert_that(maxdev!=0)
  } else {
    maxdev = 1E10 ## Some large number
  }
  assertthat::assert_that(!(is.data.frame(X)))
  assertthat::assert_that(sum(is.na(X)) == 0)
  assertthat::assert_that(length(ylist) == nrow(X))
  assertthat::assert_that(numclust > 1)
  assertthat::assert_that(niter > 1)
  ## assert_that(!(is.data.frame(ylist[[1]])))
  ## assertthat::assert_that(prob_lambda > 0)
  ## assertthat::assert_that(all(sapply(ylist, nrow) == sapply(countslist, length)))


  ## Detect if any covariates are flat (low variance)
  ## If so, remove that, run flowmix_once(), and add back zeros in the coefficients.
  flatX = which(apply(X, 2, stats::sd) < flatX_thresh)
  if(length(flatX) == 0) rm(args)
  if(length(flatX) > 0){

    warning("Some covariates are flat. The coefficients for these variables will be set to zero.")
    orig_names = colnames(X)

    ## Recursive call of flowmix_once
    args$X = X[, -flatX]
    argn <- lapply(names(args), as.name)
    names(argn) <- names(args)
    call <- as.call(c(list(as.name("flowmix_once")), argn))
    res = eval(call, args)

    if(is.null(n_princomps) & is.null(n_hidden_nodes)) {
      ## Alter the call so that those coefficients are all zero.
      res$beta <- alter_beta(res$beta, flatX, orig_names)
      res$alpha <- alter_alpha(res$alpha, flatX, orig_names)
    } else if(is.null(n_princomps)) {
      ## No PCA; direct transformation from X to X_nn, 
      ## so mark corresponding hidden layer weights as 0.
      res$W <- alter_W(res$W, flatX, orig_names)
    } 

    ## No need to set principal component factor loadings as 0 because this 
    ## is handled internally by the `prcomp` class.

    ## Alter some other things
    res$p = res$p + length(flatX)
    res$X = X

    ## Return the result
    return(res)
  }

  ## Setup
  TT = length(ylist)
  dimdat = ncol(ylist[[1]])
  p = ncol(X)

  ## Preserve original X for reporting results, since X will be 
  ## overwritting by PCA and/or neural net transformation 
  X_orig <- X

  ## Principal components transformation
  X_pc <- NULL 
  pca_obj <- NULL
  if(!is.null(n_princomps)) {
    assertthat::assert_that(is.logical(center))
    assertthat::assert_that(is.logical(scale))
    assertthat::assert_that(n_princomps %% 1 == 0)

    pca_obj <- stats::prcomp(X_orig, center = center, scale = scale)
    X_pc <- stats::predict(pca_obj, X_orig)[,seq_len(n_princomps), drop = FALSE]
    X <- X_pc
  }

  ## Neural network transformation
  X_nn <- NULL
  W <- NULL
  if(!is.null(n_hidden_nodes)) {
    assertthat::assert_that(n_hidden_nodes %% 1 == 0)

    if(length(hidden_node_dist) > 2) hidden_node_dist <- hidden_node_dist[1] # Choose "unif" by default
    assertthat::assert_that(hidden_node_dist %in% c("unif", "norm"))
    
    assertthat::assert_that(hidden_node_scale > 0)

    ## Save current seed, if it exists, or NULL if it doesn't
    prev_seed <- get0(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    ## Set user-supplied seed
    set.seed(hidden_node_seed)
    ## Randomly generate hidden layer weights
    W <- switch(hidden_node_dist, 
      unif = stats::runif(n_hidden_nodes * (ncol(X) + 1), -hidden_node_scale, hidden_node_scale), 
      norm = stats::rnorm(n_hidden_nodes * (ncol(X) + 1), 0, hidden_node_scale)
    )
    ## If RNG existed beforehand, reset seed (to prevent changing user-set RNG state)
    if(!is.null(prev_seed)) .GlobalEnv$.Random.seed <- prev_seed

    W <- matrix(W, nrow = ncol(X) + 1)
    rownames(W) <- c("intp", colnames(X))

    X_nn <- stats::plogis(cbind(1, X) %*% W)
    X <- X_nn
  }

  if(!is.null(seed)){
    assertthat::assert_that(all((seed %>% sapply(., class)) == "integer"))
    assertthat::assert_that(length(seed) == 7)
  }
  if(is.null(mn)) mn = init_mn(ylist, numclust, TT, dimdat, countslist, seed)
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
    ##if(iter == 43) browser()
    if(verbose){
      print_progress(iter-1, niter-1, "EM iterations.", start.time = start.time)
    }
    resp <- Estep(mn, sigma, prob, ylist = ylist, numclust = numclust,
                  denslist_by_clust = denslist_by_clust,
                  first_iter = (iter == 2), countslist = countslist)

    ## M step (three parts)
    ## 1. Alpha
    ## if(iter==2) browser()
    res.alpha = Mstep_alpha(resp, X, numclust, lambda = prob_lambda,
                            zerothresh = zerothresh)
    prob = res.alpha$prob
    alpha = res.alpha$alpha
    rm(res.alpha)

    ## 2. Beta
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
                                             sigma_eig_by_clust)

    ## Calculate the objectives
    objectives[iter] = objective(mn, prob, sigma, ylist,
                                 prob_lambda = prob_lambda,
                                 mean_lambda = mean_lambda,
                                 alpha = alpha, beta = beta,
                                 denslist_by_clust = denslist_by_clust,
                                 countslist = countslist)

    ## Check convergence
    converged = check_converge_rel(objectives[iter-1],
                                   objectives[iter],
                                   tol = tol_em)
    if((iter > 10) & converged){
      break()
    }
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
                        X = X_orig,
                        X_pc = X_pc,
                        n_princomps = n_princomps,
                        pca_obj = pca_obj,
                        X_nn = X_nn,
                        n_hidden_nodes = n_hidden_nodes, 
                        hidden_node_dist = hidden_node_dist, 
                        hidden_node_scale = hidden_node_scale, 
                        hidden_node_seed = hidden_node_seed,
                        W = W, 
                        prob_lambda = prob_lambda,
                        mean_lambda = mean_lambda,
                        maxdev=maxdev,
                        niter = niter,
                        admm_niters = admm_niters, 
                        seed = seed
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



##' Prediction: Given new covariates X's, generate a set of predicted cluster
##' means and probs (and return the same Sigma).
##'
##' @param object Object returned from \code{flowmix()}.
##' @param logits Logical: should the function return the logits of the cluster probs (TRUE)  
##'   or the cluster probs (FALSE)?
##' @param ... Make sure to provide new covariate values, a row vector in the
##'   same format and column names as the original X used in \code{object}, in
##'   \code{newx}.
##'
##' @return List containing mean, prob, and sigma. If 
##' user specifies logits = TRUE, prob will still be 
##' called prob, but it will actually contain the linear 
##' predictor of prob (the multinomial analogue of the 
##' logit function). 
##'
##' @export
##'
predict.flowmix <- function(object, logits = FALSE, ...){
  ## Basic checks
  ## stopifnot(ncol(new.x) == ncol(object$X))
  ## newx = X[1,,drop=FALSE]
  . = NULL ## Fixing check()
  rest = list(...)
  newx = rest[["newx"]]
  newx_pc = rest[["newx_pc"]]
  newx_nn = rest[["newx_nn"]]

  newx_list = list(newx, newx_pc, newx_nn)
  newx_is_null = sapply(newx_list, is.null)
  newx_list = newx_list[!newx_is_null]

  # Count how many newx objects are provided
  if(length(newx_list) > 1) {
    stop("Must provide only one or none of newx, newx_pc, and newx_nn.")
  } else if(length(newx_list) == 1) {
    ## Save original newx for reporting
    newx_orig = newx
    newx = newx_list[[1]]
  } else { # if all newx objects are null
    newx = object$X
    ## Save original newx for reporting
    newx_orig = newx
  }

  if(!newx_is_null[[3]]) { ## newx_nn provided
    ## Check if the variable names are the same.
    cnames = object$X_nn %>% colnames()
    cnames_new = newx %>% colnames()
    if(!all(cnames == cnames_new)) stop("Attempting to predict on X_nn but colnames do not match.")
    ## Check whether desired prediction is possible using the provided model object.
    if(is.null(object$W)) stop("Cannot predict on X_nn using a linear model.")
  } else if (!newx_is_null[[2]]) { # newx_pc provided
    ## Check if the variable names are the same.
    cnames = object$X_pc %>% colnames()
    cnames_new = newx %>% colnames()
    if(!all(cnames == cnames_new)) stop("Attempting to predict on X_pc but colnames do not match.")
    ## Check whether desired prediction is possible using the provided model object.
    if(is.null(object$pca_obj)) stop("Cannot predict on X_pc unless model has a PCA transformation.")
  } else { # newx provided or using object$X
    cnames = object$X %>% colnames()
    cnames_new = newx %>% colnames()
    ## Check whether desired prediction is possible using the provided model object.
    if(!all(cnames == cnames_new)) stop("Attempting to predict on X but colnames do not match.")
  }

  ## If newx is provided on the original X scale
  if(newx_is_null[[2]] &  newx_is_null[[3]]) {
    if(!is.null(object$pca_obj)) {
      ## If model regresses on the principal components
      ## Make the PCA transformation
      newx = stats::predict(object$pca_obj, newx)[,seq_len(object$n_princomps)]
      newx_pc = newx
    } 
  }

  ## If newx is provided on the original X scale or principal components scale
  if(newx_is_null[[3]]) {
    if(!is.null(object$W)) {
      ## If model is nonlinear
      ## Make the nonlinear transformation
      newx = stats::plogis(cbind(1, newx) %*% object$W)
      newx_nn = newx
    } 
  }

  ## Augment newx with a dummy (intercept) variable 1
  if(nrow(newx)>1){
    newx.a = cbind(rep(1, nrow(newx)), newx)
  } else {
    newx.a = c(1, newx)
  }

  TT = nrow(newx) ## This used to be nrow(X)..
  numclust = object$numclust
  dimdat = object$dimdat
  if(is.null(dimdat)) dimdat = object %>%.$mn %>% dim() %>% .[2] ## for back=compatibility

  ## Predict the means (manually).
  newmn = lapply(1:numclust, function(iclust){
    newx.a %*% object$beta[[iclust]]
  })
  newmn_array = array(NA, dim=c(TT, dimdat, numclust))
  for(iclust in 1:numclust){ newmn_array[,,iclust] = newmn[[iclust]] }

  ## Predict the probs.
  ## newprob = predict(object$alpha.fit, newx=newx, type='response')[,,1]
  probhatmat = as.matrix(tcrossprod(cbind(1, newx), object$alpha))
  if(logits) {
    newprob = probhatmat
  } else {
    probhatmat = exp(probhatmat)
    newprob = probhatmat / rowSums(probhatmat)
  }
  # probhatmat = as.matrix(exp(cbind(1,newx) %*% t(object$alpha)))
  # newprob = probhatmat / rowSums(probhatmat)

  ## predict(fit, newx=X, type="response")[,,1]
  stopifnot(all(dim(newprob) == c(TT,numclust)))
  stopifnot(logits | all(newprob >= 0)) # stop if not logits and any probs < 0 

  ## Return all three things
  return(list(mn = newmn_array,
              prob = newprob,
              pie = newprob, ## Just a copy of prob, for back-compatibility
              alpha = object$alpha,
              beta = object$beta,
              sigma = object$sigma,
              TT = object$TT,
              N = object$N,
              numclust = object$numclust,
              X = newx_orig, 
              X_pc = newx_pc, 
              X_nn = newx_nn))
}



##' Helper for making list of densities. Returns list by cluster then time
##' e.g. access by \code{denslist_by_clust[[iclust]][[tt]]}
##'
##' @inheritParams Mstep_beta_admm 
##' @param ylist T-length list each containing response matrices of size (nt x
##'   3), which contains coordinates of the 3-variate particles, organized over
##'   time (T) and with (nt) particles at every time.
##' @param mu (T x dimdat x numclust) array.
##' @param dimdat dimension of data.
##' @param numclust number of clusters.
##' @param TT number of time points
##'
##' @return numclust-lengthed list of TT-lengthed.
##'
##' @noRd
make_denslist_eigen <- function(ylist, mu,
                                TT, dimdat, numclust,
                                sigma_eig_by_clust){

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
##' @noRd
check_converge_rel <- function(old, new, tol=1E-6){ return(abs((old-new)/old) < tol )  }
