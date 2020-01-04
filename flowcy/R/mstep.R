##' Conducts the M step for estimating the coefficients for the population
##' weights.
##' @param resp Is an (T x nt x K) array.
##' @param X Covariate matrix (T x dimdat).
##' @return The multinomial logit model coefficients. A matrix of dimension (K x
##'   (p+1)).
Mstep_alpha <- function(resp, X, numclust, lambda = 0, alpha = 1,
                        refit = FALSE,
                        sel_coef = NULL,
                        bin = FALSE, ##temporary
                        thresh = 1E-8,
                        zerothresh = 1E-8
                        ){

  TT = nrow(X)
  p = ncol(X)

  ## Calculate the summed responsibilities
  resp.sum = t(sapply(resp, colSums)) ## (T x numclust)

  ## Temporary: eventually make into sparse Matrix, but for now, CVXR doesn't
  ## take sparse Matrix.
  resp.sum = as.matrix(resp.sum)
  stopifnot(dim(resp) == c(TT, numclust))

  ## Fit the model, possibly with some regularization
  if(!refit){
    ## ## TODO: wrap a try-catch around this; catch this precise message and end
    ## ## the entire algorithm by returning a flag or NULL, when this happens.
    ## fit = glmnet::glmnet(x = X, y = resp.sum, family = "multinomial",
    ##                      alpha = alpha, lambda = lambda) ## This
    ## Question: why isn't the weight of n_t becing used?
    ## ## The coefficients for the multinomial logit model are (K x (p+1))
    ## alphahat = do.call(rbind,
    ##                    lapply(coef(fit), as.numeric))

    ## Until we find an equivalence between glmnet and cvxr (unlikely), use cvxr
    ## (slow but correct):
    Xa = cbind(1, X)
    alphahat = cvxr_multinom(X = Xa, y = resp.sum, lambda = lambda, ## This was lambda=0 for no good reason
                                 exclude = 1, thresh = thresh)
    stopifnot(all(!is.na(alphahat)))
    alphahat[which(abs(alphahat) < zerothresh, arr.ind = TRUE)] = 0
    alphahat = t(as.matrix(alphahat))
    stopifnot(all(dim(alphahat) == c(numclust, (p + 1))))

  } else {
    Xa = cbind(1, X)
    alphahat = cvxr_multinom(X = Xa, y = resp.sum, lambda = lambda, ## This was lambda=0 for no good reason
                                 sel_coef = sel_coef$alpha,
                                 exclude = 1)
    alphahat[which(abs(alphahat) < zerothresh, arr.ind = TRUE)] = 0
    alphahat = t(as.matrix(alphahat))
    stopifnot(all(dim(alphahat) == c(numclust, (p + 1))))
  }


  ## While you're at it, calculate the fitted values (\pi) as well:
  piehatmat = as.matrix(exp(cbind(1, X) %*% t(alphahat)))
  piehat = piehatmat / rowSums(piehatmat)

  stopifnot(all(dim(piehat) == c(TT,numclust)))
  stopifnot(all(piehat >= 0))
  ## This should be dimension (T x K)

  return(list(pie = piehat, alpha = alphahat))
}

##' Estimates sigma. TODO: Change the format of the sigma matrix so that the
##' covariance at each time is only recorded once i.e. it is a numclust lengthed
##' list of (dimdat x dimdat) matrices.
##' @param mn are the fitted means
##' @return (TT x numclust x dimdat x dimdat) array.
Mstep_sigma_covar <- function(resp, ylist, mn, numclust,
                              bin=FALSE ){ ## Temporary

  ## Find some sizes
  TT = length(ylist)
  ntlist = sapply(ylist, nrow)
  dimdat = ncol(ylist[[1]])
  cs = c(0, cumsum(ntlist))
  irows.list = lapply(1:TT, function(tt){irows = (cs[tt] + 1):cs[tt + 1]})

  ## Set up empty residual matrix (to be reused)
  resid.rows = matrix(0, nrow = sum(ntlist), ncol=dimdat) ## trying env
  cs = c(0, cumsum(ntlist))
  vars <- vector(mode = "list", numclust)
  resp.grandsums = rowSums(sapply(1:TT, function(tt){ colSums(resp[[tt]]) }))

  for(iclust in 1:numclust){

    ## Calculate all the residuals, and weight them.
    for(tt in 1:TT){
      row.to.subtract = mn[tt,,iclust]
      resid.rows[irows.list[[tt]],] = sqrt(resp[[tt]][,iclust]) *
        sweep(ylist[[tt]], 2, rbind(row.to.subtract))
    }

    ## Calculate the covariance matrix
    myvar = crossprod(resid.rows, resid.rows) / resp.grandsums[iclust]
    vars[[iclust]] = myvar ##NEXT UP!! Work on this.
  }
  rm(resid.rows)

  ## Make into an array
  sigma_array = array(NA, dim=c(numclust, dimdat, dimdat))
  for(iclust in 1:numclust){
      sigma_array[iclust,,] = vars[[iclust]]
  }

  ## Basic check
  stopifnot(all(dim(sigma_array) == c(numclust, dimdat, dimdat)))
  return(sigma_array)
}

##' Given a matrix positive definite matrix a, computes a^{-1/2}.  Only works
##' for positive semidefinite matrices that are diagonalizable (no normal Jordan
##' forms, etc.)
##' @param a A PSD matrix.
mtsqrt_inv <- function(a){
  a.eig <- eigen(a)
  a.sqrt <- a.eig$vectors %*% diag(1 / sqrt(a.eig$values)) %*% t(a.eig$vectors)
}

##' The M step of beta, using a particular lasso regression formulation.
##' @param maxdev The desired maximum radius of the fitted means from beta0k.
##' @param sel_coef Sparsity pattern. Only active when refit=TRUE.
##' @param sigma (numclust x dimdat x dimdat) matrix.
##' @return Result of M step; a |numclust| length list of (p+1)x(d) matrices,
##'   each containing the estimated coefficients for the mean estimation.
Mstep_beta <- function(resp, ylist, X, mean_lambda=0, sigma, numclust,
                       refit = FALSE,
                       sel_coef = NULL,
                       maxdev = NULL,
                       sigma_eig_by_clust=NULL,
                       first_iter=FALSE,
                       ridge = FALSE,
                       ridge_lambda = NULL,
                       ridge_pie = NULL,
                       bin = FALSE,
                       thresh = 1E-8,
                       zerothresh = 1E-8
                       ){

  ## Preliminaries
  TT = length(ylist)
  numclust = ncol(resp[[1]])
  dimdat = ncol(ylist[[1]])
  ntlist = sapply(ylist, nrow)
  p = ncol(X)
  Xa = cbind(1, X)
  dimsigma = dim(sigma)

  if(!is.null(sigma)){
    assert_that(all.equal(dimsigma, c(numclust, dimdat, dimdat)) == TRUE)
  }

  ## Setup
  if(ridge){
    manip_obj = manip_ridge(ylist, Xa, resp, sigma, numclust,
                            sigma_eig_by_clust = sigma_eig_by_clust,
                            first_iter = first_iter,
                            ridge_lambda = ridge_lambda,
                            ridge_pie = ridge_pie)
                            ## bin = bin)
  } else {
    manip_obj = manip(ylist, Xa, resp, sigma, numclust,
                      sigma_eig_by_clust = sigma_eig_by_clust,
                      first_iter = first_iter)
                      ## bin = bin)
  }
  Xtildes = manip_obj$Xtildes
  yvecs = manip_obj$yvecs

  ## Intercepts are to be excluded from penalization.
  exclude.from.penalty = (0:(dimdat-1))*(ncol(Xa)) + 1

  ## Obtain fitted beta, separately by cluster.
  results = lapply(1:numclust, function(iclust){

    ## Give the glmnet function pre-calculated Y's and X's.
    if(is.null(maxdev) & !refit){
      ## if(refit) stop("NOT WRITTEN YET") ## Remove when done writing this
      res =  solve_lasso(x = Xtildes[[iclust]], y = yvecs[[iclust]], lambda = mean_lambda,
                         intercept = FALSE,
                         exclude.from.penalty = exclude.from.penalty)
      ## Unravel to obtain the coef and fitted response
      betahat = matrix(res$b, ncol=dimdat)
    } else {

      ## ## TEMPORARY
      ## ## First do glmnet.
      ##   res =  solve_lasso(x = Xtildes[[iclust]], y = yvecs[[iclust]], lambda = mean_lambda,
      ##                      intercept = FALSE,
      ##                      exclude.from.penalty = exclude.from.penalty)

      ## ## Unravel to obtain the coef and fitted response
      ## betahat = matrix(res$b, ncol=dimdat)
      ## slack = 1E-5 ## Just in case
      ## xb = sqrt(rowSums((X %*% betahat[-1,])^2))

      ## if(any(xb > maxdev + slack)){
      ##   print("using CVXR")

        ## If applicable, restrict the fitted means to be close to \beta0k.
      betahat = cvxr_lasso(X = Xtildes[[iclust]],
                           Xorig = X,
                           y = yvecs[[iclust]],
                           lambda = mean_lambda,
                           exclude.from.penalty = exclude.from.penalty,
                           maxdev = maxdev,
                           dimdat = dimdat,
                           numclust = numclust,
                           refit = refit,
                           sel_coef = sel_coef$beta[[iclust]],
                           thresh = thresh ## Temporary
                           )
      betahat[which(abs(betahat) < zerothresh, arr.ind = TRUE)] = 0

      ## Temporary
      ## xb = sqrt(rowSums((X %*% betahat[-1,])^2))
      ## slack = 1E-5
      ## assert_that(max(xb) <= 0.5 + slack)
      ## print(max(xb))## <= 0.5)
      ## End of Temporary

    }
    yhat = Xa %*% betahat
    assert_that(all(dim(betahat) == c(p+1,dimdat)))
    return(list(betahat = betahat, yhat = yhat))
  })

  ## Extract results; this needs to by (T x dimdat x numclust)
  betahats = lapply(results, function(a){a$betahat})
  yhats = lapply(results, function(a){as.matrix(a$yhat)})
  yhats_array = array(NA, dim = c(TT, dimdat, numclust))
  for(iclust in 1:numclust){ yhats_array[,,iclust] = yhats[[iclust]] }

  ## Each are lists of length |numclust|.
  return(list(beta = betahats,
              mns = yhats_array))
}

##' Helper to "manipulate" X and y, to get Xtilde and Ytilde and yvec for a more
##' efficient beta M step (each are |numclust|-length lists, calculated
##' separately for each cluster).
##' @param sigma only used if \code{first_iter} is TRUE. Otherwise,
##'   \code{sigma_eig_by_clust} is used.
##' @param sigma_eig_by_clust Eigendecomposition of Sigma.
##' @param first_iter TRUE if this is the first EM iteration.
##'
##' @return 3 (or dimdat) |numclust|-length lists.
manip <- function(ylist, X, resp, sigma, numclust,
                  sigma_eig_by_clust = NULL,
                  first_iter = FALSE){

  ## Make some quantities
  ntlist = sapply(ylist, nrow)
  dimdat = ncol(ylist[[1]])
  TT = nrow(X)
  resp.sum = t(sapply(resp, colSums)) ## (T x numclust)
  sigma.inv.halves = array(NA, dim=dim(sigma))

  if(first_iter){
    for(iclust in 1:numclust){
      sigma.inv.halves[iclust,,] = mtsqrt_inv(sigma[iclust,,])
    }
  } else {
    for(iclust in 1:numclust){
      sigma.inv.halves[iclust,,] = sigma_eig_by_clust[[iclust]]$inverse_sigma_half
    }
  }

  ## Pre-calculate response and covariates to feed into lasso.
  emptymat = matrix(0, nrow = dimdat, ncol = TT)
  Ytildes = list()
  for(iclust in 1:numclust){
    for(tt in 1:TT){
      emptymat[,tt] = colSums(resp[[tt]][, iclust] * ylist[[tt]])
    }
    Ytildes[[iclust]] = emptymat
  }
  Xtildes = lapply(1:numclust, function(iclust){
    sigma.inv.halves[iclust,,] %x% (sqrt(resp.sum[,iclust]) * X)
  })

  ## Vector
  yvecs = lapply(1:numclust, function(iclust){
    yvec = (1/sqrt(resp.sum[,iclust]) * t(Ytildes[[iclust]])) %*% sigma.inv.halves[iclust,,]
    yvec = as.vector(yvec)
  })

  return(list(Xtildes = Xtildes,
              Ytildes = Ytildes,
              yvecs = yvecs
              ))
}

##' Helper to "manipulate" X and y, to get Xtilde and Ytilde and yvec for a more
##' efficient beta M step (each are |numclust|-length lists, calculated
##' separately for each cluster).
##' @return 3 (or dimdat) |numclust|-length lists.
manip_ridge <- function(ylist, X, resp, sigma, numclust,
                        sigma_eig_by_clust = NULL,
                        first_iter = FALSE,
                        ridge_lambda = 0,
                        ridge_pie = NULL){
                        ## bin = FALSE
                        ## ){

  ## Make some quantities
  ntlist = sapply(ylist, nrow)
  dimdat = ncol(ylist[[1]])
  TT = nrow(X)
  ## if(bin) resp.sum = t(sapply(resp, Matrix::colSums)) ## (T x numclust)
  ## if(!bin)
  resp.sum = t(sapply(resp, colSums)) ## (T x numclust)
  sigma.inv.halves = array(NA, dim=dim(sigma))

  if(first_iter){
    for(iclust in 1:numclust){
      sigma.inv.halves[iclust,,] = mtsqrt_inv(sigma[iclust,,])
    }
  } else {
    for(iclust in 1:numclust){
      sigma.inv.halves[iclust,,] = sigma_eig_by_clust[[iclust]]$inverse_sigma_half
    }
  }

  ## Pre-calculate response and covariates to feed into lasso.
  emptymat = matrix(0, nrow = dimdat, ncol = TT)
  Ytildes = list()
  for(iclust in 1:numclust){
    for(tt in 1:TT){
      emptymat[,tt] = colSums(resp[[tt]][, iclust] * ylist[[tt]])
    }
    Ytildes[[iclust]] = emptymat
  }
  Xtildes = lapply(1:numclust, function(iclust){
    Xtilde = sigma.inv.halves[iclust,,] %x% (sqrt(resp.sum[,iclust]) * X)
    ## Append rows for the ridge
    Xridge = sqrt(ridge_lambda)  * ( diag(rep(1, dimdat)) %x% X) ##* sqrt(ridge_pie[,iclust]))
    ## if(!is.null(ridge_pie)){
    ##   Xridge =  Xridge * sqrt(ridge_pie[,iclust]) ## I think this is right.
    ## }
    return(rbind(Xtilde, Xridge))
  })

  ## Vector
  yvecs = lapply(1:numclust, function(iclust){
    yvec = (1/sqrt(resp.sum[,iclust]) * t(Ytildes[[iclust]])) %*% sigma.inv.halves[iclust,,]
    yvec = as.vector(yvec)
    ## New: augment response with zeros
    return(c(yvec, rep(0, dimdat * TT))) ## Check the dimensions
  })

  return(list(Xtildes = Xtildes,
              Ytildes = Ytildes,
              yvecs = yvecs
              ))
}


