##' The M step of beta, using ADMM. (TODO: This should be able to use the
##' eigendecomp of the Sigmas for the objective value calculation. That is next
##' up.)
##' @param niter Total number of outer iterations.
##' @param local_adapt TRUE if locally adaptive ADMM (LA-ADMM) is to be used. If
##'   so, \code{niter} becomes the inner number of iterations, and
##'   \code{local_adapt_niter} becomes the number of outer iterations.
##' @param local_adapt_niter Number of outer iterations for LA-ADMM.
##' @param space Objectives are only calculated once every \code{space} inner
##'   admm iterations.
##' @param Zs Warm-start variable.
##' @param Ws Warm-start variable.
##' @param Us Warm-start variable.
##' @param betas Warm-start variable.
##' @param rho todo: fill in
##' @param err_rel todo: fill in
##' @param err_abs todo: fill in
##' @param zerothresh todo: fill in
##'
##' @inheritParams Mstep_beta
##'
##' @return Result of M step; a |numclust| length list of (p+1)x(d) matrices,
##'   each containing the estimated coefficients for the mean estimation.
Mstep_beta_admm <- function(resp,
                            ylist,
                            X,
                            mean_lambda = 0,
                            sigma,
                            sigma_eig_by_clust = NULL,
                            first_iter = TRUE,

                            ## Warm startable variables
                            betas = NULL,
                            Zs = NULL,
                            Ws = NULL,
                            Us = NULL,
                            ## End of warm startable variables

                            maxdev = NULL,
                            niter = 1E4,
                            rho = 100, ## Some default
                            err_rel = 1E-3,
                            err_abs = 0,
                            zerothresh = 1E-6,
                            local_adapt = FALSE,
                            local_adapt_niter = 5,
                            space = 50
                            ){

  ####################
  ## Preliminaries ###
  ####################
  TT = length(ylist)
  p = ncol(X)
  numclust = ncol(resp[[1]])
  dimdat = ncol(ylist[[1]])
  ntlist = sapply(ylist, nrow)
  resp.sum = t(sapply(resp, colSums)) ## (T x numclust)
  resp.sum = as.matrix(resp.sum)
  N = sum(resp.sum) ## NEW (make more efficient, later)
  tX = t(X)
  Xa = cbind(1, X)
  Xinv = solve(t(X)%*% X + diag(1, p))
  Xaug = rbind(X, diag(1, p))

  ## Other preliminaries
  schur_syl_A_by_clust = schur_syl_B_by_clust = term3list = list()
  ybarlist = list()
  ycentered_list = Xcentered_list = yXcentered_list = list()
  Qlist = list()
  sigmainv_list = list()
  for(iclust in 1:numclust){
    ## print(iclust, numclust, "iclust")

    ## Center y and X
    obj <- weight_ylist(iclust, resp, resp.sum, ylist)
    ycentered <- obj$ycentered
    Xcentered <- center_X(iclust, resp.sum, X)
    yXcentered = ycentered %*% Xcentered
    D = diag(resp.sum[,iclust])

    ## Form the Sylvester equation coefficients in AX + XB + C = 0
    syl_A = rho * sigma[iclust,,]
    Q = 1/N * t(Xcentered) %*% D %*% Xcentered
    syl_B = Q %*% Xinv

    ## Store the Schur decomposition
    schur_syl_A_by_clust[[iclust]] = myschur(syl_A)
    schur_syl_B_by_clust[[iclust]] = myschur(syl_B)

    ## Retrieve sigma inverse from pre-computed SVD, if necessary
    if(is.null(sigma_eig_by_clust)){
      sigmainv = solve(sigma[iclust,,])
    } else {
      sigmainv = sigma_eig_by_clust[[iclust]]$sigma_inv
    }

    ## Calculate coefficients for objective value  calculation
    Qlist[[iclust]] = Q

    ## ## Also calculate some things for the objective value
    ## ylong = sweep(do.call(rbind, ylist), 2, obj$ybar)
    ## longwt = do.call(c, lapply(1:TT, function(tt){ resp[[tt]][,iclust]})) %>% sqrt()
    ## wt.long = longwt * ylong
    ## wt.ylong = longwt * ylong
    ## crossprod(wt.ylong, wt.ylong)

    ## Store the third term
    term3list[[iclust]] =  1 / N * sigmainv %*% yXcentered
    ybarlist[[iclust]] = obj$ybar

    ycentered_list[[iclust]] = ycentered
    Xcentered_list[[iclust]] = Xcentered
    yXcentered_list[[iclust]] = yXcentered
    sigmainv_list[[iclust]] = sigmainv
  }

  ##########################################
  ## Run ADMM separately on each cluster ##
  #########################################
  yhats = admm_niters = vector(length = numclust, mode = "list")
  if(first_iter) betas = vector(length = numclust, mode = "list")
  if(first_iter){
    Zs =  Ws =  Us  = vector(length = numclust, mode = "list")
  }

  fits = matrix(NA, ncol = numclust, nrow = ceiling(niter / space))

  ## For every cluster, run LA-ADMM
  resid_mat_list = list()
  start.time = Sys.time()
  for(iclust in 1:numclust){

    ## Locally adaptive ADMM.
    res = la_admm_oneclust(K = (if(local_adapt) local_adapt_niter else 1),
                           local_adapt = local_adapt,
                           iclust = iclust,
                           niter = niter,
                           p = p , TT = TT, N = N, dimdat = dimdat, maxdev = maxdev,
                           schurA = schur_syl_A_by_clust[[iclust]],
                           schurB = schur_syl_B_by_clust[[iclust]],
                           term3 = term3list[[iclust]],
                           sigmainv = sigmainv,
                           Xinv = Xinv,
                           Xaug = Xaug,
                           Xa = Xa,
                           rho = rho,
                           rhoinit = rho,
                           sigma = sigma,
                           ybar = ybarlist[[iclust]],
                           Q = Qlist[[iclust]],
                           lambda = mean_lambda,
                           resp = resp,
                           resp.sum = resp.sum,
                           ylist = ylist, X = X, tX = tX,
                           err_rel = err_rel,
                           err_abs = err_abs,
                           zerothresh = zerothresh,
                           sigma_eig_by_clust = sigma_eig_by_clust,
                           space = space,

                           ## Warm starts from previous *EM* iteration
                           first_iter = first_iter,
                           beta = betas[[iclust]],
                           U = Us[[iclust]],
                           Z = Zs[[iclust]],
                           W = Ws[[iclust]]
                           )

    ## Store the results
    betas[[iclust]] = res$beta
    yhats[[iclust]] = res$yhat
    ## fits[,iclust] = res$fits
    admm_niters[[iclust]] = res$kk

    ## Store other things for for warmstart
    Zs[[iclust]] = res$Z
    Us[[iclust]] = res$U
    Ws[[iclust]] = res$W

    ## The upper triangular matrix remains the same.
    ## The upper triangular matrix remains the same.

    resid_mat_list[[iclust]] = res$resid_mat ## temporary
  }

  ## Aggregate the yhats into one array
  yhats_array = array(NA, dim = c(TT, dimdat, numclust))
  for(iclust in 1:numclust){ yhats_array[,,iclust] = yhats[[iclust]] }

  ## Each are lists of length |numclust|.
  return(list(betas = betas,
              mns = yhats_array,
              fits = fits,
              resid_mat_list = resid_mat_list,
              admm_niters = admm_niters, ## Temporary: Seeing the number of
                                         ## outer iterations it took to
                                         ## converge.

              ## For warmstarts
              Zs = Zs,
              Ws = Ws,
              Us = Us,


              ## For using in the Sigma M step
              ycentered_list = ycentered_list,
              Xcentered_list = Xcentered_list,
              yXcentered_list = yXcentered_list,
              Qlist = Qlist
              ))
}


##' LA (locally adaptive) ADMM wrapper to \code{admm_oneclust()}.
##'
##' @param K Number of outer iterations.
##'
la_admm_oneclust <- function(K,
                             ...){

  ## Initialize arguments for ADMM.
  args <- list(...)
  p = args$p
  TT = args$TT
  dimdat = args$dimdat

  ## This initialization can come from the previous *EM* iteration.
  if(args$first_iter){
    beta = matrix(0, nrow = p + 1, ncol = dimdat)
    Z = matrix(0, nrow = TT, ncol = dimdat)
    W = matrix(0, nrow = p, ncol = dimdat)
    U = matrix(0, nrow = TT + p, ncol = dimdat)

    args[['beta']] <- beta
    args[['Z']] <- Z
    args[['W']] <- W
    args[['U']] <- U
  }

  cols = c()
  objectives = c()

  ## Temporary
  all_fits = c()
  all_cols = c()

  ## Run ADMM repeatedly with (1) double rho, and (2) previous b
  for(kk in 1:K){
    if(kk > 1){
      args[['beta']] <- beta
      args[['Z']] <- Z
      args[['W']] <- W
      args[['U']] <- U
    }

    ## Run ADMM
    if(kk > 1) args[['rho']] <- rho
    args[['outer_iter']] <- kk

    ## Call main function
    argn <- lapply(names(args), as.name)
    names(argn) <- names(args)
    call <- as.call(c(list(as.name("admm_oneclust")), argn))
    res = eval(call, args)

    ## ## Temporary (uncomment for plotting objectives)
    ## fits = as.numeric(na.omit(res$fits))
    ## all_fits = c(all_fits, fits)
    ## all_cols = c(all_cols, rep(kk, length(fits)))
    ## print(all_fits)
    ## if(length(all_fits)!=0){
    ##   plot(all_fits %>% log10(), col=all_cols, type='o', lwd=2, main = paste0("outer iter = ", kk))
    ## }
    ## ## End temporary

    objectives = c(objectives, res$fit)

    ## Handling the scenario where the objectives are all zero
    padding = 1E-12
    objectives = objectives + padding

    ## See if outer iterations should terminate
    if(outer_converge(na.omit(objectives)) | res$converge){
      break
    }

    ## Update some parameters; double the rho value
    rho = args$rho * 2
    beta = res$beta
    Z = res$Z
    W = res$W
    U = res$U

  }
  ## if(!res$converge) warning("Didn't converge at all")

  ## Record how long the admm took; in terms of # iterations.
  res$kk = kk

  return(res)
}


##' Check if 4 consecutive objective values are sufficiently close to 1.
##'
##' @param objectives Numeric vector.
##'
##' @return Boolean.
outer_converge <- function(objectives){
  consec = 4
  if(length(objectives) < consec){
    return(FALSE)
  } else {
    mytail = utils::tail(objectives, consec)
    rel_diffs = mytail[1:(consec-1)]/mytail[2:consec]


    return(all(abs(rel_diffs) - 1 < 1E-3))
  }
}



##' ADMM of one cluster.
##'
##' @param local_adapt TRUE if locally adaptive ADMM is to be used.
##'
##' @return List containing |beta|, |yhat|, |resid_mat|, |fits|.
admm_oneclust <- function(iclust, niter, Xtilde, yvec, p,
                          TT, N, dimdat, maxdev,
                          Xa,
                          rho,
                          rhoinit,
                          Xinv,
                          schurA,
                          schurB,
                          term3,
                          sigmainv,
                          Xaug,

                          ## Xa, rho,
                          ## X0, I_aug,
                          ## intercept_inds,
                          ybar,
                          Q,
                          lambda,
                          resp,
                          resp.sum,
                          ylist, X, tX, err_rel, err_abs,
                          zerothresh,
                          ## Warm startable variables
                          beta,
                          Z,
                          W,
                          U,
                          first_iter,## Not used
                          outer_iter,
                          ## End of warm startable variables
                          local_adapt,
                          sigma,
                          sigma_eig_by_clust,
                          space = 20){

  ## Initialize the variables ###
  resid_mat = matrix(NA, nrow = ceiling(niter/5), ncol = 4)
  colnames(resid_mat) = c("primresid", "primerr", "dualresid", "dualerr")
  zrows = 1:TT
  wrows = TT + (1:p)
  rhofac = rho / rhoinit

  ## Main inner LA-ADMM loop
  fits = rep(NA, ceiling(niter/space))
  converge = FALSE
  start.time = Sys.time()
  Zlist = list()

  ## This doesn't change over iterations
  schurA = myschur(schurA$orig * rhofac)
  TA = schurA$T ##* rhofac
  TB = schurB$T
  UA = schurA$Q
  UB = schurB$Q
  tUA = schurA$tQ
  tUB = schurB$tQ

  ## if(rho == 0.02) browser()
  ## print("outer_iter")
  ## print(outer_iter)
  ## print("inner iter")
  for(iter in 1:niter){
    ## print(iter)


    ## Update ytilde based on new beta
    syl_C = prepare_sylC_const3(U, Xaug, rho, Z, X, W,
                                term3,
                                sigma[iclust,,] %>% as.matrix(), Xinv)
    F = (-1) * tUA %*% syl_C %*% UB;
    beta = UA %*% matrix_function_solve_triangular_sylvester_barebones(TA, TB, F) %*% tUB
    beta = t(beta)
    if(any(is.nan(beta))) browser()

    Xbeta = X %*% beta
    ## print(summary(Xbeta))
    ## Z = Z_update(X %*% beta, U[zrows,], maxdev, rho)
    Z <- Z_updateC(Xbeta, U[zrows,, drop=FALSE], maxdev, rho, dimdat, TT)
    W = W_update(beta, U[wrows,, drop=FALSE], lambda, rho)
    if(is.vector(W)) W = cbind(W) ## make it so that it handles vectors
    U = U_update(U, rho, Xaug, beta, Z, W)
    ## print(summary(X%*%beta))

    ## Check convergence
    if( iter > 1  & iter %% 5 == 0){## & !local_adapt){

      ## Calculate convergence criterion
      obj = converge(beta, rho,
                     W, Z,
                     W_prev, Z_prev,
                     Uw = U[wrows,], Uz = U[zrows,],
                     tX = tX,
                     Xbeta = Xbeta,
                     err_rel = err_rel,
                     err_abs = err_abs)

      ## TODO change converge to use Z, W, U
      jj = (iter/ 5)
      resid_mat[jj,] = c(norm(obj$primal_resid, "F"),
                         obj$primal_err,
                         norm(obj$dual_resid,"F"),
                         obj$dual_err)

      if(obj$converge){
        converge = TRUE
        break
      }
    }

    ## ## 3. Calculate objective values for this cluster.
    W_prev = W
    Z_prev = Z
    ## Zlist[[iter]] = Z

    ## ## Temporary (uncomment for plotting objectives)
    ## ## 4. Calculate things related to convergence (Slow).
    ## space = 10
    ## if(iter %% space == 0 ){
    ##   ii = iter / space
    ##   beta0 <- intercept(resp, resp.sum, ylist, beta, X, N, iclust)
    ##   fits[ii] = objective_per_cluster(rbind(beta0, beta), ylist, Xa, resp,
    ##                                    lambda, N, dimdat, iclust, sigma, iter,
    ##                                    zerothresh, is.null(sigma_eig_by_clust),
    ##                                    sigma_eig_by_clust)
    ## }
    ## ## End of temporary
  }

  ## Gather results.
  beta = W
  beta[which(abs(beta) < zerothresh, arr.ind = TRUE)] = 0
  beta0 <- intercept(resp, resp.sum, ylist, beta, X, N, iclust,
                     ybar)

  betafull = rbind(beta0, beta)
  yhat = Xa %*% betafull

  ## if(outer_iter %% 2 == 0){

    ## ## Complete calculation (super slow)
    ## fit = objective_per_cluster(betafull, ylist, Xa, resp, lambda, N, dimdat,
    ##                             iclust, sigma, iter, zerothresh,
    ##                             is.null(sigma_eig_by_clust), sigma_eig_by_clust,
    ##                             rcpp = FALSE)


    ## ## (not needed but useful for testing) part of the objective value calculation
    ## ylong = sweep(do.call(rbind, ylist), 2, obj$ybar)
    ## longwt = do.call(c, lapply(1:TT, function(tt){ resp[[tt]][,iclust]})) %>% sqrt()
    ## wt.long = longwt * ylong
    ## wt.ylong = longwt * ylong
    ## objective_first_term = sum(diag(crossprod(wt.ylong, wt.ylong) %*% sigmainv))/ (2*N)
    ## stopifnot(fit, objective_first_term + sum(diag(Q)) - sum(diag(M)) + lambda * sum(abs(beta) > zerothresh))

    ## ## Todo: get rid of the transposes
    Q2 = (1 / 2) * Q %*% (beta %*% sigmainv %*% t(beta))
    M = t(term3) %*% t(beta)
    fit = sum(diag(Q2)) - sum(diag(M)) + lambda * sum(abs(beta) > zerothresh)

  ## } else {
  ##   fit = NA
  ## }

  return(list(beta = betafull,
              yhat = yhat,
              resid_mat = resid_mat,
              fits = fits,
              converge = converge,
              fit = fit,
              ## Other variables to return.
              Z = Z,
              W = W,
              U = U
              ))
}


##' Update U in the new ADMM
U_update <- function(U, rho, Xaug, beta, Z, W){
  (U + rho * (Xaug %*% beta - rbind(Z, W)) )
}

##' Update U in the new ADMM
W_update  <- function(beta, Uw, lambda, rho){
  soft_thresh(beta + Uw/rho, lambda/rho)
}

##' (Helper) Soft thresholding of |a| at radius of |b|.
##'
##' @param a Numeric vector.
##' @param b Numeric vector of same length as \code{length(b)}.
##'
##' @return Soft-thresholded version of a and b.
##'
soft_thresh <- function(a, b){
  return(sign(a) * pmax(0, abs(a)-b))

  ## ## (Kind of) test for this function:
  ## x = seq(from=-10, to=10, length=100)
  ## plot(y=soft_thresh(x, 3), x=x)
  ## abline(h=0)
  ## abline(v=c(-3,3), col='grey')
  ## abline(v=0)
}


intercept <- function(resp, resp.sum, ylist, beta, X, N, iclust, ybar){
  dimdat = ncol(ylist[[1]])
  TT = length(ylist)
  resp.sum.thisclust = sum(resp.sum[,iclust])
  wt.resid.sum = rep(0, dimdat)
  mn = X %*% beta
  yhat = (resp.sum[,iclust] / resp.sum.thisclust) * mn
  return(ybar - colSums(yhat))
}


Z_update  <- function(Xbeta, Uz, C, rho){
  mat = Xbeta + Uz/rho
  Z = projCmat(mat, C)
}

projCmat <- function(mat, C){
  if(!is.null(C)){
    vlens = sqrt(rowSums(mat * mat))
    inds = which(vlens > C)
    if(length(inds) > 0){
      mat[inds,] = mat[inds,] * C / vlens[inds]
    }
  }
  return(mat)
}


###### All helpers for main ADMM function, prepping for sylvester solver ###


center_X <- function(iclust, resp.sum, X){
  resp.sum.thisclust = sum(resp.sum[,iclust])
  Xtilde = colSums(resp.sum[,iclust] * X) / resp.sum.thisclust
  Xcentered = sweep(X, 2, Xtilde, check.margin=FALSE)
  return(Xcentered)
}


##' Not used anymore; erase soon
center_y <- function(iclust, ylist, resp, resp.sum){
  resp.sum.thisclust = sum(resp.sum[,iclust])
  ybar = Reduce("+", Map(function(y, myresp){
    myresp[,iclust, drop = TRUE] %*% y
  }, ylist, resp)) / resp.sum.thisclust
  ycentered_list = Map(function(y, myresp){
    colSums(myresp[,iclust] * sweep(y, 2, ybar, check.margin=FALSE))
  }, ylist, resp) ##/ resp.sum.thisclust
  ycentered = do.call(cbind, ycentered_list)
  return(ycentered)
}

weight_ylist <- function(iclust, resp, resp.sum, ylist){

  ## Setup
  dimdat = ncol(ylist[[1]])
  TT = length(ylist)

  ## All weighted data
  weighted_ylist = Map(function(myresp, y){
    myresp[,iclust, drop = TRUE] * y
  }, resp, ylist)

  ## All weighted data SUMs
  weighted_ysum = lapply(weighted_ylist, colSums)

  ## Grand mean of data
  resp.sum.thisclust = sum(resp.sum[,iclust])
  ybar = Reduce("+", weighted_ysum) / resp.sum.thisclust

  ## Centered weighted ylist
  weighted_ybar = resp.sum[,iclust] * matrix(ybar, TT, dimdat, byrow=TRUE)
  ycentered = do.call(rbind, weighted_ysum) - weighted_ybar##sweep(do.call(cbind, centered_y), 1, ybar)
  return(list(weighted_ylist = weighted_ylist,
              weighted_ysum = weighted_ysum,
              ybar = ybar,
              ycentered = t(ycentered)))
}

##' Convenience function for obtaining a Schur decomposition of a matrix
##' \code{mat}. Used for simplifying a Sylvester equation.
##'
##' @param mat Numeric matrix.
##'
##' @return List.
myschur <- function(mat){
  stopifnot(nrow(mat) == ncol(mat))
  if(is.numeric(mat) & length(mat)==1) mat = mat %>% as.matrix()
  obj = Matrix::Schur(mat)
  obj$tQ = t(obj$Q)
  obj$orig = mat
  return(obj)
}



##' (Helper) Check convergence for ADMM.
##'
##' @param err_rel = 1E-3
##' @param inheritParams Mstep_beta_admm
##'
##'
converge <- function(beta1, rho, w, Z, w_prev, Z_prev, Uw, Uz, tX, Xbeta1,
                     err_rel = 1E-4,
                     err_abs = 0
                     ){

  ## Form primal and dual residuals, and other things.
  prim1 = rbind(beta1, Xbeta1) ## (I don't like the rbind)
  prim2 = - rbind(w, Z)
  primal_resid = prim1 + prim2
  dual_resid = -rho * ((w - w_prev) + (tX %*% (Z - Z_prev)))
  tAU = Uw + tX %*% Uz

  ## Form primal and dual tolerances.
  primal_err = sqrt(length(primal_resid)) * err_abs +
    err_rel * max(norm(prim1, "F"), norm(prim2, "F"))
  dual_err = sqrt(length(dual_resid)) * err_abs +
    err_rel * norm(tAU, "F")

  ## Check convergence.
  primal_resid_size = norm(primal_resid, "F")
  dual_resid_size = norm(dual_resid, "F")
  primal_converge = ( primal_resid_size  <= primal_err )
  dual_converge = ( dual_resid_size <= dual_err )

  ## Some checks (trying to address problems with |converge|).
  assertthat::assert_that(is.numeric(primal_resid_size))
  assertthat::assert_that(is.numeric(primal_err))
  assertthat::assert_that(is.numeric(dual_resid_size))
  assertthat::assert_that(is.numeric(dual_err))

  ## return(primal_converge & dual_converge)
  converge = primal_converge & dual_converge
  if(is.na(converge)) browser()

  return(list(primal_resid = primal_resid,
              primal_err = primal_err,
              dual_resid = dual_resid,
              dual_err = dual_err,
              converge = converge))
}


##' (Helper) calculates the per-cluster objective value for the ADMM. Super
##' inefficient right now.
##'
##' @param beta p x d matrix
##' @inheritParams Mstep_beta_admm
##'
objective_per_cluster <- function(beta, ylist, Xa, resp, lambda, N, dimdat,
                                  iclust, sigma, iter, zerothresh, first_iter,
                                  sigma_eig_by_clust=NULL, rcpp = FALSE){

  ## Setup
  TT = length(ylist)
  p = ncol(Xa) - 1
  mn = Xa %*% beta  ## VERY SLOW (600)

  ## Obtain square root of sigma.
  if(first_iter){
    sigma_half <- sigma_half_from_eig(eigendecomp_sigma_barebones(sigma[iclust,,]))
  } else {
    sigma_eig <- sigma_eig_by_clust[[iclust]]
    sigma_half <- sigma_eig$inverse_sigma_half
  }


  ## Obtain the weighted residuals once.
  ## save(ylist, TT, mn, resp, iclust, sigma_half, file="~/Desktop/rcpp-speed.Rdata")
  ## load(file="~/Desktop/rcpp-speed.Rdata", verbose=TRUE)
  if(!rcpp){
      wt_resids_list = lapply(1:TT, function(tt){
        y = ylist[[tt]]
        mumat = matrix(mn[tt,],
                       ncol = ncol(y),
                       nrow = nrow(y),
                       byrow = TRUE)
        wt_resid = sqrt(resp[[tt]][, iclust]) * (y - mumat)
      })
      transformed_resids = do.call(rbind, wt_resids_list) %*% sigma_half
      resid.prods = rowSums(transformed_resids * transformed_resids)
  }
  if(rcpp){


    ## This works and is slightly faster
    ## wt_resids_list = list()
    ## for(tt in 1:TT){
    ##   wt_resids_list[[tt]] = subtractC3(sqrt(resp[[tt]][, iclust]),
    ##                                     ylist[[tt]],
    ##                                     mn[tt,])
    ## }
    form_wt_resid <- function(tt){
      subtractC3(sqrt(resp[[tt]][, iclust]),
                 ylist[[tt]],
                 mn[tt,])
    }
    wt_resids_list = lapply(1:TT, form_wt_resid)
    transformed_resids = do.call(rbind, wt_resids_list) %*% sigma_half
    resid.prods = rowSumsC2_arma(transformed_resids)

    ## Restructuring of data.
    ## ylong = do.call(rbind, ylist) ## This can be done in \code{flowmix()}, once.
    ## ntlist = sapply(ylist, nrow)
    ## mumat = mn[rep(1:TT, ntlist), ] ## You can't get around this. ## this takes a long time.
    ## longwt = do.call(c, lapply(1:TT, function(tt){ resp[[tt]][,iclust]})) %>% sqrt()
    ## resid.prods = dothisC(longwt, ylong, mumat, sigma_half)
  }

  ## ## Calculate the inner product: resid * sigma^-1 * t(resid)
  grand.sum = sum(resid.prods)

  ## Calculate the objective value
  obj = (1/(2 * N)) * grand.sum + lambda * sum(abs(beta[-1,]) > zerothresh) ## Another bug!
  return(obj)
}
