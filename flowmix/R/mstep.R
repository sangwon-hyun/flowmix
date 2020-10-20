##' Solves the M step for estimating the coefficients for the cluster
##' probabilities.
##'
##' @param resp Responsibilities; an (T x nt x K) array.
##' @param X Covariate matrix (T x dimdat).
##'
##' @return The multinomial logit model coefficients. A matrix of dimension (K x
##'   (p+1)).
Mstep_alpha <- function(resp, X, numclust, lambda,
                        zerothresh = 1E-8){

  ## Basic checks
  TT = nrow(X)
  p = ncol(X)

  ## Calculate the summed responsibilities
  resp.sum = t(sapply(resp, colSums)) ## (T x numclust)
  resp.sum = as.matrix(resp.sum)
  stopifnot(dim(resp) == c(TT, numclust))

  ## Fit the model
  alpha = NULL

  ## Try glmnet first:
  tryCatch({
    alpha = alpha.glmnet = solve_multinom(resp.sum, X, lambda)
  }, error=function(e){})

  ## In the rare case where the above fails, use CVXR:
  if(is.null(alpha)){
    alpha = alpha.cvxr = cvxr_multinom(resp.sum, cbind(1,X), lambda,
                                       exclude.from.penalty = 1,
                                       N = sum(resp.sum))
  }
  ## Two notes
  ## 1. A unit test exists to compare the likelihood between \code{glmnet} andd
  ## \code{CVXR} solutions. Here is a simple one, commented out for now:
  ## objective.cvxr = multinom_objective(alpha.cvxr, cbind(1,X), resp.sum, lambda, N=sum(resp.sum))
  ## objective.glmnet = multinom_objective(alpha.glmnet, cbind(1,X), resp.sum, lambda, N=sum(resp.sum))
  ## stopifnot(abs(objective.cvxr - objective.glmnet) < 1E-3)
  ## 2. The estimated alphas are not unique, only unique up to a constant shift.
  ## See \url{https://www.jstatsoft.org/article/view/v033i01} section 4.1 for more details.

  ## Threshold some of the alpha values.
  alpha[which(abs(alpha) < zerothresh, arr.ind = TRUE)] = 0
  alpha = t(as.matrix(alpha))
  stopifnot(all(dim(alpha) == c(numclust, (p + 1))))

  ## Calculate the fitted values (\pi) as well:
  Xa = cbind(1, X)
  piehatmat = as.matrix(exp(Xa %*% t(alpha)))
  piehat = piehatmat / rowSums(piehatmat)

  ## ## Start of temporary

  ## alpha2[which(abs(alpha2) < zerothresh, arr.ind = TRUE)] = 0
  ## alpha2 = t(as.matrix(alpha2))
  ## stopifnot(all(dim(alpha2) == c(numclust, (p + 1))))

  ## ## Calculate the fitted values (\pi) as well:
  ## Xa = cbind(1, X)
  ## piehatmat = as.matrix(exp(Xa %*% t(alpha2)))
  ## piehat2 = piehatmat / rowSums(piehatmat)

  ## print(range(piehat2 - piehat))
  ## ## End of temporary


  ## Checking dimensions one last time.
  stopifnot(all(dim(piehat) == c(TT,numclust)))
  if(any(is.na(piehat))){ print(piehat); stop("piehat was erroneous") }
  if(!(all(piehat >=0))){ print(piehat); stop("piehat was erroneous") }
  stopifnot(all(piehat >= 0))

  return(list(pie = piehat, alpha = alpha))
}


##' Estimate sigma.
##'
##' @param mn Fitted means.
##' @param resp Responsibilities.
##' @param ylist List of cytograms.
##' @param numclust Number of clusters.
##'
##' @return An array of size (numclust x dimdat x dimdat) containing the
##'   covariance matrices.
Mstep_sigma_covar_new <- function(resp, ylist, mn, numclust,
                                  beta,
                                  ycentered_list,
                                  Xcentered_list,
                                  yXcentered_list,
                                  Qlist
                                  ){
                                  ## resp, ylist, mn, numclust){

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
  ylong = do.call(rbind, ylist)
  ntlist = sapply(ylist, nrow)
  ## resp.grandsums = rowSums(sapply(1:TT, function(tt){ colSums(resp[[tt]]) }))

  for(iclust in 1:numclust){

    ## Calculate all the residuals, and weight them.
    ## Xa = cbind(1,X)
    ## mn2 = Xa %*% beta[[iclust]]##[-1,]
    ## mn[,,iclust] - mn2

    ## Calculate things the hard way
    resp.thisclust = lapply(resp, function(myresp) myresp[,iclust, drop = TRUE])
    resp.long = do.call(c, resp.thisclust)
    irows = rep(1:nrow(mn), times = ntlist)
    mnlong = mn[irows,,iclust]
    resid.rows = sqrt(resp.long) * (ylong - mnlong)

    ## Also try it the "alternate" way with centered y's and X's.
    resp.sum = t(sapply(resp, colSums)) ## (T x numclust)
    resp.sum = as.matrix(resp.sum)
    resp.thisclust = lapply(resp, function(myresp) myresp[,iclust, drop = TRUE])
    resp.long = do.call(c, resp.thisclust)
    irows = rep(1:nrow(mn), times = ntlist)
    mnlong = mn[irows,,iclust]
    resids = (ylong - mnlong)
    obj <- weight_ylist(iclust, resp, resp.sum, ylist)
    ycentered_long = sweep(ylong, 2, obj$ybar)
    Xcentered <- center_X(iclust, resp.sum, X)
    this_beta = beta[[iclust]][-1,]
    mn_from_centered = Xcentered %*% this_beta
    resid.rows2 = sqrt(resp.long) * (ycentered_long - mn_from_centered[rep(1:TT, times=ntlist),])
    stopifnot(all(abs(resid.rows - resid.rows2)<1E-10))


    ## Calculate the covariance matrix
    myvar = crossprod(resid.rows, resid.rows) / sum(resp.long)##resp.grandsums[iclust]

    ## save(ylong, mn, resp, file="~/Desktop/temp-mstep.Rdata")
    ## load( file="~/Desktop/temp-mstep.Rdata")

    ## Now, first, load quantities that've been formed already in beta M step
    N = sum(unlist(resp)) ## NEW (make more efficient, later)
    this_beta = beta[[iclust]][-1,]
    yXcentered = yXcentered_list[[iclust]]
    ycentered = ycentered_list[[iclust]]
    Xcentered = Xcentered_list[[iclust]]
    Q = Qlist[[iclust]]

    ## (Instead of the above blocK: manually create them instead)
    resp.sum = t(sapply(resp, colSums)) ## (T x numclust)
    resp.sum = as.matrix(resp.sum)
    obj <- weight_ylist(iclust, resp, resp.sum, ylist)
    ## ycentered <- center_y(iclust, ylist, resp, resp.sum)
    ycentered <- obj$ycentered
    Xcentered <- center_X(iclust, resp.sum, X)
    yXcentered = ycentered %*% Xcentered
    D = diag(resp.sum[,iclust])
    C = sum(unlist(resp.thisclust))

    ## 1. I've checked the correctness of this
    resp.thisclust = lapply(resp, function(myresp) myresp[,iclust, drop = TRUE])
    resp.long = do.call(c, resp.thisclust)
    ylong = sqrt(resp.long) * do.call(rbind, ylist)
    first_term = crossprod(ylong, ylong)
    ## resid.rows = matrix(0, nrow = sum(ntlist), ncol=dimdat) ## trying env
    ## for(tt in 1:TT){ resid.rows[irows.list[[tt]],] = sqrt(resp[[tt]][,iclust]) * ylist[[tt]] }
    ## myvar0 = crossprod(resid.rows, resid.rows) / resp.grandsums[iclust]
    ## stopifnot(all(myvar0 == first_term/C))


    ## 2. How to check the correctness of this?
    ## second_term0 = (yXcentered %*% this_beta)
    ## second_term = second_term0 + t(second_term0)

    ## Checking ycentered and yXentered
    weighted_ylist_centered = Map(function(myresp, y){
      myresp[,iclust, drop = TRUE] * sweep(y,2, obj$ybar)
    }, resp, ylist)
    weighted_ysum_centered = lapply(weighted_ylist_centered, colSums)
    ycentered2 = do.call(rbind, weighted_ysum_centered)
    stopifnot(all.equal(t(ycentered2), obj$ycentered))
    ycentered = obj$ycentered
    yXcentered2 = ycentered %*% Xcentered
    stopifnot(all.equal(yXcentered2, yXcentered))

    second_term0 = (yXcentered %*% this_beta)
    second_term = second_term0 + t(second_term0)

    ## 3. How to check the correctness of this?
    ## third_term = t(this_beta) %*% (Q * N) %*% this_beta
    third_term = t(this_beta) %*% t(Xcentered) %*% D %*% Xcentered %*% this_beta
    mnhat_centered  = Xcentered %*% this_beta

    ## Putting it all together
    myvar_new = (first_term - second_term + third_term) / C
    ## first_term/C
    ## second_term/C
    ## third_term/C

    ## Testing it
    myvar_new/myvar
    myvar %>% round(3)
    myvar_new %>% round(3)



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



##' Estimate sigma.
##'
##' @param mn Fitted means.
##' @param resp Responsibilities.
##' @param ylist List of cytograms.
##' @param numclust Number of clusters.
##'
##' @return An array of size (numclust x dimdat x dimdat) containing the
##'   covariance matrices.
Mstep_sigma_covar <- function(resp, ylist, mn, numclust){

  ## Find some sizes
  TT = length(ylist)
  ntlist = sapply(ylist, nrow)
  dimdat = ncol(ylist[[1]])
  cs = c(0, cumsum(ntlist))
  irows.list = lapply(1:TT, function(tt){irows = (cs[tt] + 1):cs[tt + 1]})

  ## Set up empty residual matrix (to be reused)
  cs = c(0, cumsum(ntlist))
  vars <- vector(mode = "list", numclust)
  ylong = do.call(rbind, ylist)
  ntlist = sapply(ylist, nrow)
  irows = rep(1:nrow(mn), times = ntlist)

  for(iclust in 1:numclust){
      resp.thisclust = lapply(resp, function(myresp) myresp[,iclust, drop = TRUE])
      resp.long = do.call(c, resp.thisclust)
      mnlong = mn[irows,,iclust]
      vars[[iclust]] = estepC(ylong, mnlong, sqrt(resp.long), sum(resp.long))
  }

  ## Make into an array
  sigma_array = array(NA, dim=c(numclust, dimdat, dimdat))
  for(iclust in 1:numclust){
      sigma_array[iclust,,] = vars[[iclust]]
  }

  ## Basic check
  stopifnot(all(dim(sigma_array) == c(numclust, dimdat, dimdat)))
  return(sigma_array)
}

##' Given a matrix positive definite matrix a, compute a^{-1/2}.  Only works for
##' positive semidefinite matrices that are diagonalizable (no normal Jordan
##' forms, etc.)
##'
##' @param a A PSD matrix.
##'
##' @return Matrix of the same size as \code{a}.
mtsqrt_inv <- function(a){
  a.eig <- eigen(a)

  ## In case vec is a single element, in which case diag() isn't quite right.
  vec = 1 / sqrt(a.eig$values)
  if(length(vec)==1){
    mat = vec
  } else {
    mat = diag(vec)
  }

  ## a.sqrt <- a.eig$vectors %*% diag(1 / sqrt(a.eig$values)) %*% t(a.eig$vectors)
  a.sqrt <- a.eig$vectors %*% mat %*% t(a.eig$vectors)
}

##' Solves the M step of beta, using a particular lasso regression
##' formulation. The preferred default over this function is
##' \code{Mstep_admm()}, which is much faster.
##'
##' @param maxdev The desired maximum radius of the fitted means from beta0k.
##' @param sel_coef Sparsity pattern. Only active when refit=TRUE. List
##'   containing boolean matrices named "alpha" and "beta".
##' @param sigma (numclust x dimdat x dimdat) matrix.
##'
##' @return Result of M step; a |numclust| length list of (p+1)x(d) matrices,
##'   each containing the estimated coefficients for the mean estimation.
Mstep_beta <- function(resp, ylist, X, mean_lambda = 0, sigma,
                       refit = FALSE,
                       sel_coef = NULL,
                       maxdev = NULL,
                       sigma_eig_by_clust = NULL,
                       first_iter = FALSE,
                       ridge = FALSE,
                       ridge_lambda = NULL,
                       ridge_pie = NULL,
                       cvxr_ecos_thresh = 1E-8,
                       cvxr_scs_eps = 1E-5,
                       zerothresh = 1E-8
                       ){

  ## Preliminaries
  TT = length(ylist)
  numclust = ncol(resp[[1]])
  dimdat = ncol(ylist[[1]])
  ntlist = sapply(ylist, nrow)
  ## N = sum(ntlist) ## BUG!!
  resp.sum = t(sapply(resp, colSums)) ## (T x numclust)
  resp.sum = as.matrix(resp.sum)
  N = sum(resp.sum) ## NEW (make more efficient, later)
  p = ncol(X)
  Xa = cbind(1, X)
  dimsigma = dim(sigma)

  if(!is.null(sigma)){
    assertthat::assert_that(all.equal(dimsigma, c(numclust, dimdat, dimdat)) == TRUE)
  }

  ## Setup
  if(ridge){
    manip_obj = manip_ridge(ylist, Xa, resp, sigma, numclust,
                            sigma_eig_by_clust = sigma_eig_by_clust,
                            first_iter = first_iter,
                            ridge_lambda = ridge_lambda,
                            ridge_pie = ridge_pie)
  } else {
    manip_obj = manip(ylist, Xa, resp, sigma, numclust,
                      sigma_eig_by_clust = sigma_eig_by_clust,
                      first_iter = first_iter)
  }
  Xtildes = manip_obj$Xtildes
  yvecs = manip_obj$yvecs

  ## Intercepts are to be excluded from penalization.
  exclude.from.penalty = (0:(dimdat-1)) * (ncol(Xa)) + 1

  ## Obtain fitted beta, separately by cluster.
  results = lapply(1:numclust, function(iclust){

    betahat = cvxr_lasso(X = Xtildes[[iclust]],
                         Xorig = X,
                         y = yvecs[[iclust]],
                         lambda = mean_lambda,
                         exclude.from.penalty = exclude.from.penalty,
                         maxdev = maxdev,
                         dimdat = dimdat,
                         N = N,
                         refit = refit,
                         sel_coef = sel_coef$beta[[iclust]],
                         ecos_thresh = cvxr_ecos_thresh,
                         scs_eps = cvxr_scs_eps)## temporary
    betahat[which(abs(betahat) < zerothresh, arr.ind = TRUE)] = 0

    ## ## Double checking
    ## xb = sqrt(rowSums((X %*% betahat[-1,])^2))
    ## slack = 1E-2
    ## assert_that(max(xb) <= 0.5 + slack)

    yhat = Xa %*% betahat
    assertthat::assert_that(all(dim(betahat) == c(p+1,dimdat)))
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
##'
##' @return 3 (or dimdat) |numclust|-length lists.
manip_ridge <- function(ylist, X, resp, sigma, numclust,
                        sigma_eig_by_clust = NULL,
                        first_iter = FALSE,
                        ridge_lambda = 0,
                        ridge_pie = NULL){

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

  Xunderbar = diag(rep(1, dimdat)) %x% X
  Xtildes = lapply(1:numclust, function(iclust){
    Xtilde = sigma.inv.halves[iclust,,] %x% (sqrt(resp.sum[,iclust]) * X)
    ## Append rows for the ridge
    Xridge = sqrt(ridge_lambda)  * Xunderbar
    return(rbind(Xtilde, Xridge))
  })

  ## Vector
  yvecs = lapply(1:numclust, function(iclust){
    yvec = (1/sqrt(resp.sum[,iclust]) * t(Ytildes[[iclust]])) %*% sigma.inv.halves[iclust,,]
    yvec = as.vector(yvec)
    ## New: augment response with zeros
    return(c(yvec, rep(0, dimdat * TT)))
  })

  return(list(Xtildes = Xtildes,
              Ytildes = Ytildes,
              yvecs = yvecs
              ))
}






##' Helper function to calculate multinomial objective. This calculates:
##'
##'    \eqn{ \hat \alpha \leftarrow \argmax_{\alpha_{0k}, \alpha_k}
##'      \frac{1}{N}\sum_{t=1}^T \left( \sum_{k=1}^K \gamma_{\cdot kt}
##'      (\alpha_{0k} + {X^{(t)}}^T \alpha_k) - n_t \log \sum_{l=1}^K
##'      \exp(\alpha_{0l} + {X^{(t)}}^T \alpha_l) \right) - \lambda_\alpha
##'      \sum_{k=1}^K \|\alpha_k\|_1}
##'
##' @param alpha (p x dimdat) matrix of the alpha coefficients.
##' @param x (T x p) covariate matrix.
##' @param y (T x dimdat) matrix.
##' @param lambda regularization parameter value.
##' @param N scaling factor for the likelihood function.
##' @param exclude.from.penality If not NULL, contains index in 1:p
##'
##' @return Penalized objective value.
multinom_objective <- function(alpha, x, y, lambda, N,
                               exclude.from.penalty=NULL) {
  n <- nrow(x)
  p <- ncol(x)
  L <- ncol(y)
  eta <- x %*% alpha
  v = 1:p
  if(!is.null(exclude.from.penalty)){
    stopifnot(all(exclude.from.penalty %in% (1:p)))
    v = (1:p)[-exclude.from.penalty]
 }
  ys = rowSums(y)
  (1/N) * (sum(eta * y) - sum(ys * log(rowSums(exp(eta))))) -
    lambda * sum(abs(alpha[v,]))
}
