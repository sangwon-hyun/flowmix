## Synopsis: these are the specific, scale-consistent lasso solvers used
## internally in this package.

##' Solving the  no-intercept Lasso problem using  CVXR.  \deqn{\min_{\beta} 1/2n
##' \|y - X\beta\|^2 + \lambda \|\beta\|_1}
##' @param X Covariate matrix.
##' @param y Response vector.
##' @param lambda regularization problem.
cvxr_lasso <- function(y, X, lambda, exclude.from.penalty=NULL, thresh=1E-12){
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
  obj <- sum_squares(y - X %*% beta) / (2 * n) + lambda * sum(abs(beta[v]))
  prob <- Problem(Minimize(obj))
  result <- solve(prob, FEASTOL = thresh, RELTOL = thresh, ABSTOL = thresh)
  return(as.numeric(result$getValue(beta)))
}

##' This solves the lasso problem: \deqn{\min_{\beta_0, \beta} 1/2\|y-1\beta_0-
##' X\beta\|^2 + \lambda \|\beta\|_1}. This is the interface function; all other
##' helpers are not meant to be used by end-user.
##' @param y Response vector (n vector).
##' @param X Covariate matrix (n by p). Not assumed to be standardized.
##' @param intercept TRUE if intercept is desired. Defaults to TRUE.
##' @return list that contains linear coefficients (length p), intercept
##'   (scalar, or NULL if intercept=TRUE), and fitted values (length n).
solve_lasso <- function(y, x, lambda, intercept=TRUE, exclude.from.penalty=NULL){

  ## Setup
  n = nrow(x)
  p = ncol(x)

  ## Solve lasso problem with or without intercept.
  if(intercept){
    if(!is.null(exclude.from.penalty)){
      stop("Doesn't support exclude.from.penalty option when intercept exists!")
    }
    res = glmnet(y=y, x=x, lambda=lambda, intercept=TRUE, standardize=FALSE)
    b = as.numeric(coef(res))[2:(p+1)]
    b0 = as.numeric(coef(res))[1]
  } else {
    pen = rep(1, p)
    if(!is.null(exclude.from.penalty)){ pen[exclude.from.penalty] = 0
      pen = pen / sum(pen) * (p)
      lam = lambda / unique(pen[which(pen!=0)])
    }
    res = glmnet(y=y, x=x, lambda=lam, intercept=FALSE, standardize=FALSE,
                 penalty.factor=pen)
    b0 = NULL
    b = as.numeric(coef(res))[2:(p+1)]
  }
  yhat = predict(res, newx=x)
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
