##' The M step of beta, using ADMM. (TODO: This should be able to use the
##' eigendecomp of the Sigmas for the objective value calculation. That is next
##' up.)
##'
##' @param maxdev The desired maximum radius of the fitted means from beta0k.
##' @param sel_coef Sparsity pattern. Only active when refit=TRUE.
##' @param sigma (numclust x dimdat x dimdat) matrix.
##' @param local_adapt TRUE if locally adaptive ADMM (LA-ADMM) is to be used. If
##'   so, \code{niter} becomes the inner number of iterations, and
##'   \code{local_adapt_niter} becomes the number of outer iterations.
##' @param local_adapt_niter Number of outer iterations for LA-ADMM.
##'
##' @return Result of M step; a |numclust| length list of (p+1)x(d) matrices,
##'   each containing the estimated coefficients for the mean estimation.
Mstep_beta_admm_new <- function(resp,
                            ylist,
                            X,
                            mean_lambda = 0,
                            sigma,
                            sigma_eig_by_clust = NULL,
                            first_iter = TRUE,
                            ## em_iter,

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
                            plot.admm = FALSE,
                            local_adapt = FALSE,
                            local_adapt_niter = 5,
                            space = 50,
                            rcpp = FALSE
                            ){

  ####################
  ## Preliminaries ###
  ####################
  TT = length(ylist)
  p = ncol(X)
  numclust = ncol(resp[[1]])
  dimdat = ncol(ylist[[1]])
  ntlist = sapply(ylist, nrow)
  ## N = sum(ntlist) ## OLD BUG!!
  resp.sum = t(sapply(resp, colSums)) ## (T x numclust)
  resp.sum = as.matrix(resp.sum)
  N = sum(resp.sum) ## NEW (make more efficient, later)
  tX = t(X)
  Xa = cbind(1, X)
  Xinv = solve(t(X)%*% X + diag(1, p))
  Xaug = rbind(X, diag(1, p))

  ## Other preliminaries
  schur_syl_A_by_clust = schur_syl_B_by_clust = term3list = list()
  for(iclust in 1:numclust){
    ## print(iclust, numclust, "iclust")

    ## Center y and X
    ycentered <- center_y(iclust, ylist, resp, resp.sum)
    Xcentered <- center_X(iclust, resp.sum, X)
    yXcentered = ycentered %*% Xcentered
    D = diag(resp.sum[,iclust])

    ## Form the Sylvester equation coefficients in AX + XB + C = 0
    syl_A = rho * sigma[iclust,,]
    syl_B = 1/N * t(Xcentered) %*% D %*% Xcentered %*% Xinv

    ## Store the Schur decomposition
    schur_syl_A_by_clust[[iclust]] = myschur(syl_A)
    schur_syl_B_by_clust[[iclust]] = myschur(syl_B)

    ## Retrieve sigma inverse from pre-computed SVD, if necessary
    if(is.null(sigma_eig_by_clust)){
      sigmainv = solve(sigma[iclust,,])
    } else {
      sigmainv = sigma_eig_by_clust[[iclust]]$sigma_inv
    }

    ## Store the third term
    term3list[[iclust]] = 1 / N * sigmainv %*% yXcentered
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

    ## temporary LA-ADMM.
    res = la_admm_oneclust_new(K = (if(local_adapt) local_adapt_niter else 1),
                           local_adapt = local_adapt,
                           iclust = iclust,
                           niter = niter,
                           p = p , TT = TT, N = N, dimdat = dimdat, maxdev = maxdev,

                           ## Also important: the Schur factorization of A and B in AX + XB + C = 0
                           schurA = schur_syl_A_by_clust[[iclust]],
                           schurB = schur_syl_B_by_clust[[iclust]],
                           term3 = term3list[[iclust]],
                           ## sigmainv = sigmainv,
                           sigmainv = solve(sigma[iclust,,]),
                           Xinv = Xinv,

                           Xaug = Xaug,
                           ## Xa = Xa,
                           ## rho = rho,
                           ## X0 = X0, I_aug = I_aug,
                           Xa = Xa,
                           rho = rho,
                           rhoinit = rho,
                           sigma = sigma,
                           lambda = mean_lambda,
                           resp = resp,
                           resp.sum = resp.sum,
                           ylist = ylist, X = X, tX = tX,
                           err_rel = err_rel,
                           err_abs = err_abs,
                           zerothresh = zerothresh,
                           sigma_eig_by_clust = sigma_eig_by_clust,
                           plot = plot.admm,
                           space = space,
                           rcpp = rcpp,

                           ## Warm starts from previous *EM* iteration
                           first_iter = first_iter,
                           ## em_iter = em_iter,
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
              Us = Us
              ))
}


##' LA (locally adaptive) ADMM wrapper to \code{admm_oneclust()}.
##'
##' @param K Number of outer iterations.
##'
la_admm_oneclust_new <- function(K,
                             ...){

  ## Initialize arguments for ADMM.
  args <- list(...)
  p = args$p
  TT = args$TT
  dimdat = args$dimdat

  ## This initialization can come from the previous *EM* iteration.
  if(args$first_iter){
    beta = matrix(0, nrow=p+1, ncol=dimdat)
    ## Z = matrix(0, nrow = TT, ncol = dimdat)
    ## wvec = rep(0, p * dimdat)
    ## uw  = rep(0, p * dimdat)
    ## Uz = matrix(0, nrow = TT, ncol = dimdat)

    ## New
    Z = matrix(0, nrow = TT, ncol = dimdat)
    W = matrix(0, nrow = p, ncol = dimdat)
    U = matrix(0, nrow = TT + p, ncol = dimdat)
    ## End of new


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
    call <- as.call(c(list(as.name("admm_oneclust_new")), argn))
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


    ## See if outer iterations should terminate
    objectives = c(objectives, res$fit)
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
outer_converge_new <- function(objectives){
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
admm_oneclust_new <- function(iclust, niter, Xtilde, yvec, p,
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
                          ## em_iter,
                          ## End of warm startable variables
                          local_adapt,
                          sigma,
                          sigma_eig_by_clust,
                          space = 20,
                          plot = FALSE,
                          rcpp = FALSE){

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
  print("outer_iter")
  print(outer_iter)
  print("inner iter")
  for(iter in 1:niter){
    print(iter)
    ## printprogress(iter, niter, "iters")
    ## printprogress(iter, "iters")

    ## Update ytilde based on new beta
    syl_C = prepare_sylC_const3(U, Xaug, rho, Z, X, W,
                                term3,
                                sigma[iclust,,], Xinv)
    F = (-1) * tUA %*% syl_C %*% UB;
    beta = UA %*% matrix_function_solve_triangular_sylvester_barebones(TA, TB, F) %*% tUB
    beta = t(beta)
    if(any(is.nan(beta))) browser()

    ## ##### Beginning of Test version ################3
    ## ycentered <- center_y(iclust, ylist, resp, resp.sum)
    ## Xcentered <- center_X(iclust, resp.sum, X)
    ## yXcentered = ycentered %*% Xcentered
    ## D = diag(resp.sum[,iclust])

    ## ## Form the Sylvester equation coefficients in AX + XB + C = 0
    ## syl_B = 1/N * t(Xcentered) %*% D %*% Xcentered %*% Xinv
    ## syl_A = rho * sigma[iclust,,]

    ## ## Retrieve sigma inverse from pre-computed SVD, if necessary
    ## sigmainv = solve(sigma[iclust,,])

    ## ## Store the third term
    ## term3 = 1 / N * sigmainv %*% yXcentered
    ## syl_C = prepare_sylC_const3(U, Xaug, rho, Z, X, W,
    ##                             term3,
    ##                             sigma[iclust,,], Xinv)
    ## beta = sylC(syl_A, syl_B, syl_C) %>% t()


    ## syl_A_reconstructed = UA %*% TA %*% tUA
    ## syl_B_reconstructed = UB %*% TB %*% tUB
    ## syl_A - syl_A_reconstructed
    ## dim(syl_A)
    ## syl_A_reconstructed
    ## sylC(




    ## browser()
    ## beta0 <- intercept(resp, resp.sum, ylist, betanew, X, N, iclust)
    ##### End of Test version ########################


    Xbeta = X %*% beta
    ## print(summary(Xbeta))
    ## Z = Z_update(X %*% beta, U[zrows,], maxdev, rho)
    Z <- Z_updateC(Xbeta, U[zrows,], maxdev, rho, dimdat, TT)
    W = W_update(beta, U[wrows,], lambda, rho)
    U = U_update(U, rho, Xaug, beta, Z, W)
    ## print(summary(X%*%beta))



    ## Check convergence
    if( iter > 1  & iter %% 5 == 0){## & !local_adapt){

      ## Calculate convergence criterion
      obj = converge(beta, rho,
                     W, Z,
                     w_prev, Z_prev,
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
        print("converged!")
        converge = TRUE
        break
      }
    }

    ## ## 3. Calculate objective values for this cluster.
    w_prev = W
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
  beta[which(abs(beta) < zerothresh, arr.ind = TRUE)] = 0
  beta0 <- intercept(resp, resp.sum, ylist, beta, X, N, iclust)
  betafull = rbind(beta0, beta)
  yhat = Xa %*% betafull

  if(outer_iter %% 2 == 0){
    fit = objective_per_cluster(betafull, ylist, Xa, resp, lambda, N, dimdat,
                                iclust, sigma, iter, zerothresh,
                                is.null(sigma_eig_by_clust), sigma_eig_by_clust,
                                rcpp = rcpp)
  } else {
    fit = NA
  }

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
              ## wvec = wvec,
              ## uw = uw,
              ## Uz = Uz
              ))
}

##' Update beta using sylvester equations.
beta_update <- function(beta, ylist, rho, sigma,
                        X, Xinv, Xaug,
                        iclust,
                        U, Z, W,
                        yXcentered,
                        ## Xcentered_Xinv,
                        syl_B,
                        D,
                        N,
                        sigmainv
                        ){

  ## Setup
  dimdat = ncol(ylist[[1]])

  ## Solve Sylvester Equation. (SH: correct)
  syl_A = rho * sigma[iclust,,]
  ## syl_B = 1/N * t(Xcentered) %*% D %*% Xcentered_Xinv ##Xcentered %*% Xinv
  ## syl_C = sigma[iclust,,] %*% C2 %*% Xinv
  syl_C = prepare_sylC_const3(U, Xaug, rho, Z, X, W, N, sigmainv, yXcentered, sigma[iclust,,], Xinv)



  sol = t(sylC(syl_A, syl_B, syl_C))## %>% t()

  return(sol)
}


##' Update U in the new ADMM
U_update <- function(U, rho, Xaug, beta, Z, W){
  U + rho * (Xaug %*% beta - rbind(Z, W))
}

##' Update U in the new ADMM
W_update  <- function(beta, Uw, lambda, rho){
  soft_thresh(beta + Uw/rho, lambda/rho)
}

intercept <- function(resp, resp.sum, ylist, beta, X, N, iclust){
  dimdat = ncol(ylist[[1]])
  TT = length(ylist)
  wt.resid.sum = rep(0, dimdat)
  mn = X %*% beta

  ## ## Is this any faster if we do a complete restructuring?
  ## resp_long = do.call(c,  lapply(1:TT, function(tt){resp[[tt]][,iclust, drop=TRUE]}))
  ## ylong = do.call(rbind, ylist)
  ## mnlong = do.call(rbind, mnlong)

  for(tt in 1:TT){
    resid_tt = sweep(ylist[[tt]], 2,
                     mn[tt,,drop=TRUE], check.margin=FALSE)
    resp_tt = resp[[tt]][,iclust, drop=TRUE]
    wt.resid = resp_tt * resid_tt
    wt.resid.sum = wt.resid.sum + colSums(wt.resid)
  }

  ## TODO There must be a faster way to do this.
  return(wt.resid.sum / sum(resp.sum[,iclust]))##N)
}



sigma_half_from_eig_temp <- function(sigma_eig){
  vec = sqrt(sigma_eig$values)
  if(length(vec)==1){
    mat = vec
  } else {
    mat = diag(vec)
  }
  (sigma_eig$vectors %*% mat %*% t(sigma_eig$vectors))
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
  ## basic check (commented out for now)
  ## stopifnot(resp.sum.this.clust == sum(resp.sum[,iclust]))
  resp.sum.thisclust = sum(resp.sum[,iclust])
  Xtilde = colSums(resp.sum[,iclust] * X) / resp.sum.thisclust
  Xcentered = sweep(X, 2, Xtilde, check.margin=FALSE)
  return(Xcentered)
}
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


myschur <- function(mat){
  stopifnot(nrow(mat) == ncol(mat))
  obj = Matrix::Schur(mat)
  obj$tQ = t(obj$Q)
  obj$orig = mat
  return(obj)
}




