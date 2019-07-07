## Synopsis: these are the specific, scale-consistent lasso solvers used
## internally in this package.

##' Solving the  no-intercept Lasso problem using  CVXR.  \deqn{\min_{\beta} 1/2n
##' \|y - X\beta\|^2 + \lambda \|\beta\|_1}
##' @param X Covariate matrix.
##' @param y Response vector.
##' @param lambda regularization problem.
cvxr_lasso <- function(y, X, lambda, exclude.from.penalty=NULL, thresh=1E-12
                       ## coef_limit=NULL
                       ){
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

  ## if(!is.null(coef_limit)){## Experimental feature
  ##   constraints <- list(beta[v] >= min(coef_limit),
  ##                       beta[v] <= max(coef_limit))
  ## } else {
    constraints = list()
  ## }

  ## Perform elastic-net regression
  obj <- sum_squares(y - X %*% beta) / (2 * n) + lambda * sum(abs(beta[v]))
  prob <- Problem(Minimize(obj), constraints)
  result <- solve(prob, FEASTOL = thresh, RELTOL = thresh, ABSTOL = thresh)
  return(as.numeric(result$getValue(beta)))
}

##' This solves the lasso problem: \deqn{\min_{\beta_0, \beta} 1/2\|y-1\beta_0-
##' X\beta\|^2 + \lambda \|\beta\|_1}. There is an option
##' \code{exclude.from.penalty} to remove some coordinates of \eqn{beta} from
##' penalization.
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
    res = glmnet::glmnet(y=y, x=x, lambda=lambda, intercept=TRUE, standardize=FALSE)
    b = as.numeric(coef(res))[2:(p+1)]
    b0 = as.numeric(coef(res))[1]
  } else {
    pen = rep(1, p)
    lam = lambda
    if(!is.null(exclude.from.penalty)){
      pen[exclude.from.penalty] = 0
      pen = pen / sum(pen) * (p)
      lam = lambda / unique(pen[which(pen!=0)])
    }
    res = glmnet::glmnet(y=y, x=x, lambda=lam, intercept=FALSE, standardize=FALSE,
                         penalty.factor=pen
                         ## lower.limits=min(coef_limit),## Experimental feature
                         ## upper.limits=max(coef_limit)## Experimental feature
                         )
    b0 = NULL
    b = as.numeric(coef(res))[2:(p+1)]
  }
  yhat = predict(res, newx=x)
  return(list(b=b, b0=b0, yhat=yhat))
}
