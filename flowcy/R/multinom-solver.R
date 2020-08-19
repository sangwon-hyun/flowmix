## Synopsis: these are the specific, scale-consistent sparse multinom solvers
## used internally in this package.


##' Solves the l1-penalized multinom problem using \code{glmnet} (not assuming
##' that $y$ have row sums of 1):
##'
##' \deqn{\frac{1}{n} \sum_{i=1}^n \sum_{k=1}^K y_{ik} ( \alpha_{0k} +
##' \alpha_k^T X^{(t)}) - \log \left( \sum_{l=1}^K y_{ik} \exp ( \alpha_{0l} +
##' \alpha_l^T X^{(t)}) \right) - \lambda \sum_{l=1}^K \left(\| \alpha_l \|_1
##' \right)}
##'
##' Note, this should have exactly the same output as \code{cvxr_multinom()}
##' (using \code{cbind(1, X)} instead of \code{X}, and
##' \code{exclude.from.penalty=1}); \code{cvxr_multinom()} is intended as a
##' fall-back strategy for when glmnet fails, and as an internal tester.
##'
##' @param y response (TT by numclust).
##' @param X Covariates (TT by p).
##' @param lambda regularization parameter for l1 penalization.
##' @param lambda_max Default is \code{10}. Internally, \code{solve_multinom()}
##'   uses glmnet on a logarithmically spaced decreasing sequence of lambda
##'   values, whose last entry is \code{lambda}.
##'
##' @return A (p+1) by (numclust) matrix.
solve_multinom <- function(y, X, lambda, lambda_max = 10){

  ## Basic check
  if(lambda > lambda_max){
    lambda_max = lambda * 100
  }

  TT = nrow(X)
  ysums = rowSums(y)
  N = sum(ysums)
  stopifnot(lambda > 0)
  lambdas = exp(seq(from = log(10), to = log(lambda), length = 30))
  fit <- glmnet::glmnet(x = X,
                        y = y/ysums,
                        lambda = lambdas,
                        family = "multinomial",
                        intercept = TRUE,
                        weights = ysums / N * TT)
  coefs = glmnet::coef.glmnet(fit, s = lambda)
  return(as.matrix(do.call(cbind, coefs)))
}


##' Solves, using CVXR, the l1-penalized multinom problem (not assuming that $y$
##' have row sums of 1)
##'
##'      \deqn{\hat \alpha \leftarrow \argmax_{, \alpha_k}
##'      \frac{1}{N}\sum_{t=1}^T \left( \sum_{k=1}^K \gamma_{\cdot kt}
##'      ({X^{(t)}}^T \alpha_k) - n_t \log \sum_{l=1}^K
##'      \exp({X^{(t)}}^T \alpha_l) \right) - \lambda_\alpha
##'      \sum_{k=1}^K \|\alpha_k[-v]\|_1}.
##'
##' Here, \deqn{v \in \{1, \dots, p\}} are indices that should be excluded from
##' the l1 penalty, specified in \code{exclude_from_penalty}. One example of its
##' usage is for \code{X} to be \code{cbind(1,X)}, and
##' \code{exclude_from_penalty=1}, to consider a model with an intercept, but
##' not penalizing the intercept coefficient.
##'
##' @param y Matrix valued response; each row is an observation, and each column
##'   is a discrete outcome out of K (n x K).
##' @param X Covariate matrix (n x p)
##' @param lambda Regularization parameter for l1 penalized estimation.
##' @param thresh ECOS solver threshold.
##' @param N Total number of particles in all cytograms i.e. \deqn{N=\sum_t n_t}.
##'
##' @return (p x K) coefficient matrix.
cvxr_multinom <- function(y, X, lambda,
                          exclude.from.penalty=NULL, thresh = 1E-8,
                          N){

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
  obj1 = CVXR::sum_entries(CVXR::mul_elemwise(y, (X %*% alphamat))) ## Using version 99.7 or before
  ## obj1 = CVXR::sum_entries(CVXR::multiply(y, (X %*% alphamat))) ## Using version 1.0 or above; this needs some work, since the
  ##                                                               ## \code{solve} gives an S4 related error that is so far unresolved.

  ## Second component
  ry = rowSums(y)
  obj2 <- sum((CVXR::log_sum_exp( X %*% alphamat, 1)) * ry)

  ## Sum them
  ## obj <- (obj1  - obj2)  / TT - lambda * sum(abs(alphamat[v,]))
  obj <- (obj1  - obj2)  / N - lambda * sum(abs(alphamat[v,])) ## Putting things in the scale of a single particle.

  ## Solve the problem using the default, ECOS CVXR solver.
  prob <- CVXR::Problem(CVXR::Maximize(obj))
  result <- solve(prob, solver="ECOS",
                  FEASTOL = thresh, RELTOL = thresh, ABSTOL = thresh)

  ## If all goes well, return the optimizer.
  alphamat = result$getValue(alphamat)
  alphamat[which(abs(alphamat) < 1E-8)] = 0
  return(alphamat)
}
