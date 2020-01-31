## Synopsis: these are the specific, scale-consistent sparse multinom solvers
## used internally in this package.

##' Solves the l1-penalized multinom problem using GLMNET (not assuming that $y$
##' have row sums of 1):
##'
##' \deqn{\frac{1}{n} \sum_{i=1}^n \sum_{k=1}^K y_{ik} ( \alpha_{0k} +
##' \alpha_k^T X^{(t)}) - \log \left( \sum_{l=1}^K y_{ik} \exp ( \alpha_{0l} +
##' \alpha_l^T X^{(t)}) \right) - \lambda \sum_{l=1}^K \left(\| \alpha_l \|_1
##' \right)}
##'
##' Note, this should have exactly the same output as cvxr_multinom() (using
##' cbind(1, X) instead of X, and exclude.from.penalty=1) which is an internal
##' tester.
##'
##' @param y response (TT by numclust).
##' @param X Covariates (TT by p).
##' @param lambda regularization parameter for l1 penalization.
##'
##' @return A (p+1) by (numclust) matrix.
solve_multinom <- function(y, X, lambda){

  TT = nrow(X)
  ysums = rowSums(y)
  N = sum(ysums)
  fit <- glmnet::glmnet(x = X,
                        y = y/ysums,
                        lambda = lambda,
                        family = "multinomial",
                        intercept = TRUE,
                        weights = ysums / N * TT)
  return(as.matrix(do.call(cbind,coef(fit))))
}

##' Solves, using CVXR, the l1-penalized multinom problem (not assuming that $y$
##' have row sums of 1):
##'
##' \deqn{\frac{1}{n} \sum_{i=1}^n \sum_{k=1}^K y_{ik} ( \alpha_{0k} +
##' \alpha_k^T X^{(t)}) - \log \left( \sum_{l=1}^K y_{ik} \exp ( \alpha_{0l} +
##' \alpha_l^T X^{(t)}) \right) - \lambda \sum_{l=1}^K \left(\| \alpha_l \|_1
##' \right)}
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
  obj1 = CVXR::sum_entries(CVXR::mul_elemwise(y, (X %*% alphamat)))

  ## Second component
  obj2 <- sum((CVXR::log_sum_exp( X %*% alphamat, 1)) * rowSums(y))

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
