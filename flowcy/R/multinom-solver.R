## Synopsis: these are the specific, scale-consistent multinom lasso solvers used
## internally in this package.

##' (DOESN'T solve right problem) Solves the l1-penalized multinom problem:
##'
##' \deqn{ \frac{1}{n} \sum_{i=1}^n
##' \sum_{k=1}^K y_{ik} (\alpha_{0k} + x_i^T \alpha_k)
##' - (\sum{k=1^K} y_{ik}) \log \left( \sum_{l=1}^K \exp ( \alpha_{0l} + x^T \alpha_l) \right)
##' - \lambda \sum_{l=1}^K \left(\| \alpha_l \|_1 \right)}
##'
##' When there is no intercept, \eqn{\alpha_{0k}} is just set to zero.
##' @param y response (TT by numclust).
##' @param X Covariates (TT by p).
##' @param lambda regularization parameter for l1 penalization.
##' @param intercept TRUE if intercept should be included in model. Intercept is
##'   not regularized.
##' @return A (p+1) by (numclust) matrix.
solve_multinom <- function(y, X, lambda,
                           coef_limit=NULL, ## Experimental feature
                           alpha=1, ##  temporary
                           ntlist=NULL
                           ){
  if(is.null(ntlist))ntlist=rep(1,nrow(y))

  ysums = colSums(y)
    fit <- glmnet::glmnet(x = X,
                          y = y/ysums,
                          lambda = TT / sum(wts) * lambda,
                          family="multinomial",
                          intercept=TRUE,
                          weights=ntlist ## temporary
                          ## alpha=alpha
                          ## lower.limits=min(coef_limit),## Experimental feature
                          ## upper.limits=max(coef_limit)## Experimental feature
                          ## exclude=exclude ## Experimental feature
                          )
    return(as.matrix(do.call(cbind,coef(fit))))
}

##' Solves, using CVXR, the l1-penalized
##' multinom problem:
##'
##' \deqn{\frac{1}{n} \sum_{i=1}^n \sum_{k=1}^K y_{ik} (
##' x_i^T \alpha_k) - \log \left( \sum_{l=1}^K \exp ( x^T \alpha_l) \right) -
##' \lambda \sum_{l=1}^K \left(\| \alpha_l \|_1 \right)}
##'
##' @param y Matrix valued response; each row is an observation, and each column
##'   is a discrete outcome out of K (n x K).
##' @param X Covariate matrix (n x p)
##' @param lambda Regularization parameter for l1 penalized estimation.
##' @param thresh ECOS solver threshold.
##'
##' @return (p x K) coefficient matrix.
cvxr_multinom <- function(y, X, lambda, exclude.from.penalty=NULL, thresh = 1E-8){

  ## Setup
  numclust = ncol(y)
  TT = nrow(X)
  p = ncol(X)
  v = 1:p
  if(!is.null(exclude.from.penalty)){
    stopifnot(all(exclude.from.penalty %in% (1:p)))
    v = (1:p)[-exclude.from.penalty]
  }
  alphamat <- CVXR::Variable(p, numclust)

  ## First component
  obj1 = CVXR::sum_entries(CVXR::mul_elemwise(y, (X %*% alphamat)))

  ## Second component
  obj2 <- sum((CVXR::log_sum_exp( X %*% alphamat, 1)) * rowSums(y))

  ## Sum them
  obj <- (obj1  - obj2)  / TT - lambda * sum(abs(alphamat[v,]))

  ## Solve the problem using one of two CVXR solvers.
  prob <- CVXR::Problem(CVXR::Maximize(obj))
  result <- solve(prob, solver="ECOS",
                  FEASTOL = thresh, RELTOL = thresh, ABSTOL = thresh)

  ## If all goes well, return the optimizer.
  return(result$getValue(alphamat))
}
