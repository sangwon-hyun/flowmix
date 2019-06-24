## Synopsis: these are the specific, scale-consisten lasso solvers used in this
## package.

##' Solving the  no-intercept Lasso problem using  CVXR.  \deqn{\min_{\beta} 1/2n
##' \|y - X\beta\|^2 + \lambda \|\beta\|_1}
##' @param X Covariate matrix.
##' @param y Response vector.
##' @param lambda regularization problem.
cvx_lasso <- function(y, X, lambda, exclude.from.penalty=NULL, thresh=1E-12){
  n = nrow(X)
  p = ncol(X)
  beta <- Variable(p)
  loss <- sum((y - X %*% beta)^2) / (2 * n)

  ## Set up exclusion from penalty, if applicable.
  if(is.null(exclude.from.penalty)){
    v = 1:p
  } else {
    assert_that(all(exclude.from.penalty %in% (1:p)))
    v = (1:p)[-exclude.from.penalty]
  }

  ## Perform elastic-net regression
  obj <- sum_squares(y - X %*% beta) / (2 * n) + lambda * p_norm(beta[v], 1)
  prob <- Problem(Minimize(obj))
  result <- solve(prob, FEASTOL = thresh, RELTOL = thresh, ABSTOL = thresh)
  as.numeric(result$getValue(beta))
}


##' Same as cvx_lasso().  Results should match exactly.  Solving the
##' no-intercept Lasso problem using CVXR.  \deqn{\min_{\beta} 1/2n \|y -
##' X\beta\|^2 + \lambda \|\beta\|_1}
glmnet_lasso <- function(y, X, lambda, exclude.from.penalty=NULL){
  lambdascaled = sqrt(sum(y*y)) * lambda
  Xscaled = X * sqrt(sum(y*y))
  penalty.factor = rep(1, ncol(X))
  assert_that(all(exclude.from.penalty %in% (1:p)))
  if(!is.null(exclude.from.penalty)) penalty.factor[exclude.from.penalty] = 0
  fit = glmnet::glmnet(x=Xscaled, y=y, lambda=lambdascaled, alpha=1,
                       intercept=FALSE, standardize=FALSE,
                       penalty.factor=penalty.factor)
  b = coef(fit)[-1,]
  b = b * sqrt(sum(y*y))
  return(as.numeric(b))
}


##' This solves the lasso  problem: \deqn{\min_{\beta_0, \beta} 1/2\|y-1\beta_0-
##' X\beta\|^2 + \lambda \|\beta\|_1}. This is the interface function; all other
##' helpers are not meant to be used by end-user.
##' @param y Response vector.
##' @param X Covariate matrix. Assumed to be standardized.
##' @param intercept TRUE if intercept is desired. Defaults to TRUE.
##' @return list that contains coefficients, fitted values, and glmnet object.
solve_lasso <- function(y, X, lambda, intercept=TRUE, exclude.from.penalty=NULL){

  ## Setup
  n = nrow(X)
  p = ncol(X)

  ## Basic Checks
  mncheck = apply(X, 2, function(mycol){
    mean(mycol)==0
  })
  sdcheck = apply(X, 2, function(mycol){
    sd(mycol)==1
  })

  if(intercept){

    ## Solve the Lasso problem with intercept.
    Xcentered = X - colMeans(X)
    ycentered = y - mean(y)
    ## b = glmnet_lasso(ycentered, Xcentered, lambda, exclude.from.penalty)
    b = cvx_lasso(Xcentered, ycentered, lambda, exclude.from.penalty)
    b0 = mean(y) - colMeans(X) %*% b
    yhat = rep(b0, n) + X %*% b

  } else {

    ## Solve the Lasso problem /without/ intercept.
    ## b = glmnet_lasso(y, X, lambda, exclude.from.penalty)
    b = cvx_lasso(y, X, lambda, exclude.from.penalty)
    b0 = NULL
    yhat = X %*% b

  }
  return(list(b=b, b0=b0, yhat=yhat))
}


##' Solves the l1-penalized multinom problem: ##' \deqn{ \frac{1}{n} \sum_{i=1}^n \sum_{k=1}^K y_{ik} (\alpha_{0k} + x_i^T
##' \alpha_k) -
##' \log \left( \sum_{l=1}^K \exp ( \alpha_{0l} + x^T \alpha_l) \right)
##' - \lambda \sum_{l=1}^K \left(\| \alpha_l \|_1 \right)}
##' When there is no intercept, \eqn{\alpha_{0k}} is just set to zero.
##' @param y response (TT by numclust).
##' @param X Covariates (TT by p).
##' @param lambda regularization parameter for l1 penalization.
##' @param intercept TRUE if intercept should be included in model. Intercept is
##'   not regularized.
##' @return A (p+1) by (numclust) matrix.
solve_multinom <- function(y, X, lambda, intercept){
    fit <- glmnet(x=X,
                  y=y,
                  lambda=lambda,
                  family="multinomial",
                  intercept=intercept)
    return(do.call(cbind,coef(fit)))
}

##' Solves, using CVXR, the l1-penalized multinom problem:
##' \deqn{ \frac{1}{n} \sum_{i=1}^n \sum_{k=1}^K y_{ik} ( x_i^T \alpha_k) -
##' \log \left( \sum_{l=1}^K \exp ( x^T \alpha_l) \right)
##' - \lambda \sum_{l=1}^K \left(\| \alpha_l \|_1 \right)}
cvxr_multinom <- function(y, X, lambda, exclude.from.penalty=NULL){
  ## Setup
  numclust = ncol(y)
  p = ncol(X)
  v = 1:p
  if(!is.null(exclude.from.penalty)){
    stopifnot(all(exclude.from.penalty %in% (1:p)))
    v = (1:p)[-exclude.from.penalty]
  }
  alphamat <- Variable(p, numclust)

  ## First component
  obj1 = matrix_trace(t(y) %*% X %*% alphamat) / TT## X

  ## Second component
  obj2 <- sum(log_sum_exp(X %*% alphamat, 1)) / TT

  ## Sum them
  obj <- obj1  - obj2 - lambda * sum(abs(alphamat[v,]))

  ## Solve the problem
  prob <- Problem(Maximize(obj))
  result <- solve(prob)
  return(result$getValue(alphamat))
}
