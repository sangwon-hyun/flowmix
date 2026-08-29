##' Solves the M step for estimating the coefficients for the cluster
##' probabilities.
##'
##' Note, the estimated alphas are not unique, only unique up to a constant
##' shift.  See \url{https://www.jstatsoft.org/article/view/v033i01} section 4.1
##' for more details.
##'
##' @param resp Responsibilities; an (T x nt x K) array.
##' @param X Covariate matrix (T x dimdat).
##' @param lambda Regularization parameter.
##' @param zerothresh Values below \code{zerothresh} are set to zero.
##' @param numclust Number of clusters.
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
  alpha = alpha.glmnet = solve_multinom(resp.sum, X, lambda)

  ## Threshold some of the alpha values.
  alpha[which(abs(alpha) < zerothresh, arr.ind = TRUE)] = 0
  alpha = t(as.matrix(alpha))
  stopifnot(all(dim(alpha) == c(numclust, (p + 1))))

  ## Calculate the fitted values (\pi) as well:
  Xa = cbind(1, X)
  probhatmat = as.matrix(exp(Xa %*% t(alpha)))
  probhat = probhatmat / rowSums(probhatmat)

  ## Checking dimensions one last time.
  stopifnot(all(dim(probhat) == c(TT,numclust)))
  if(any(is.na(probhat))){ print(probhat); stop("probhat was erroneous") }
  if(!(all(probhat >=0))){ print(probhat); stop("probhat was erroneous") }
  stopifnot(all(probhat >= 0))

  return(list(prob = probhat, alpha = alpha))
}




##' M-step for covariance matrix \eqn{\Sigma_k} for cluster \eqn{k} in 1 through
##' \code{numclust}.
##'
##' @param mn Fitted means.
##' @param resp Responsibilities.
##' @param ylist List of cytograms.
##' @param numclust Number of clusters.
##'
##' @return An array of size (numclust x dimdat x dimdat) containing the
##'   covariance matrices.
Mstep_sigma <- function(resp, ylist, mn, numclust){

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
      if(is.vector(mnlong)) mnlong = mnlong %>% cbind()
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

##' Given a matrix positive definite matrix a, compute \eqn{a^{-1/2}}.  Only works for
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
##' @noRd
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
