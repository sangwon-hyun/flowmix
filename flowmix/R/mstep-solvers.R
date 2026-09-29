##' Solves the l1-penalized multinom problem using \code{glmnet} (not assuming
##' that $y$ have row sums of 1):
##'
##' \deqn{\frac{1}{n} \sum_{i=1}^n \sum_{k=1}^K y_{ik} ( \alpha_{0k} +
##' \alpha_k^T X^{(t)}) - \log \left( \sum_{l=1}^K y_{ik} \exp ( \alpha_{0l} +
##' \alpha_l^T X^{(t)}) \right) - \lambda \sum_{l=1}^K \left(\| \alpha_l \|_1
##' \right)}
##'
##' @param y response (TT by numclust).
##' @param X Covariates (TT by p).
##' @param lambda regularization parameter for l1 penalization.
##' @param lambda_max Default is \code{10}. Internally, \code{solve_multinom()}
##'   uses glmnet on a logarithmically spaced decreasing sequence of lambda
##'   values, whose last entry is \code{lambda}.
##'
##' @return A (p+1) by (numclust) matrix.
##' @noRd
solve_multinom <- function(y, X, lambda, lambda_max = 10){

  ## Basic check
  if(lambda > lambda_max){
    lambda_max = lambda * 100
  }

  TT = nrow(X)
  ysums = rowSums(y)
  N = sum(ysums)
  stopifnot(lambda > 0)
  lambdas = exp(seq(from = log(lambda_max), to = log(lambda), length = 30))

  fit = tryCatch({
    glmnet::glmnet(x = X,
                          y = y/ysums,
                          lambda = lambdas,
                          family = "multinomial",
                          intercept = TRUE,
                          weights = ysums / N * TT)
  }, error = function(e){return(NULL)})

  ## If an error is thrown (e.g., happens when ys are equal over all time
  ## points) try adding slight amount of noise to the counts.
  if(is.null(fit)){
    numclust = ncol(y)
    y = y + stats::rnorm(length(y), 0, 0.01)
    ysums = rowSums(y)
    fit = glmnet::glmnet(x = X,
                          y = y/ysums,
                          lambda = lambdas,
                          family = "multinomial",
                          intercept = TRUE,
                          weights = ysums / N * TT)
  }
  coefs = glmnet::coef.glmnet(fit, s = lambda)
  return(as.matrix(do.call(cbind, coefs)))
}

