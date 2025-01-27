##' Solve the no-intercept lasso problem using CVXR. This is not a general lasso
##' function, but designed for the beta M step of the EM algorithm. The response
##' is a (dimdat T) vector, and the covariate matrix is a ((dimdat T) x (dimdat
##' p)). This solves the optimization and organizes+returns the optimization
##' variable as a ((p+1) x dimdat) matrix.
##'
##' \deqn{\min_{\beta} 1/2 \|y - X\vect(\beta)\|^2 + \lambda \|\beta\|_1}
##'
##' In addition, this includes other modifications such as the ball constraint
##' on the means over time, whose radius is limited to be \code{maxdev}.
##'
##' Internally, it tries an ECOS solver first, then it tries a SCS solver. The
##' latter is supposedly faster (and less exact) for larger problems, but I
##' haven't seen that in our application FYI.
##'
##' @param X Covariate matrix.
##' @param y Response vector.
##' @param lambda Regularization parameter.
##' @param maxdev Radius of ball constraint. Defaults to NULL.
##'
##' @return Fitted beta matrix.
##' @noRd
##'
##' @importFrom CVXR Variable sum_squares vec square Minimize solve
cvxr_lasso <- function(y, X,  lambda, Xorig=NULL,
                       exclude.from.penalty = NULL,
                       thresh = 1E-8,
                       maxdev = NULL,
                       dimdat,
                       N,
                       ecos_thresh = 1E-8,
                       scs_eps = 1E-5){

  ## Define dimensions
  n = nrow(X)
  p = ncol(X)
  TT = (n/dimdat)
  pp = p/dimdat - 1 ## The minus 1 is for the intercept

  ## Define the parameter
  betamat <- CVXR::Variable(rows = pp+1,
                            cols = dimdat)

  ## Set up exclusion from penalty, if applicable.
  if(is.null(exclude.from.penalty)){
    v = 1:p
  } else {
    assertthat::assert_that(all(exclude.from.penalty %in% (1:p)))
    v = (1:p)[-exclude.from.penalty]
  }

  ## Set up l1-penalized regression.
  obj = CVXR::sum_squares(y - X %*% CVXR::vec(betamat)) / (2 * N)
  obj = obj + lambda * sum(abs(CVXR::vec(betamat)[v]))

  ## Setup the Xbeta penalty.
  constraints = list()
  if(!is.null(maxdev)){
    constraints = list(CVXR::sum_entries(CVXR::square(Xorig %*% betamat[-1,]), 1) <= rep(maxdev^2, TT))
  }

  ## Try all two CVXR solvers.
  prob <- CVXR::Problem(CVXR::Minimize(obj), constraints)
  result = NULL
  result <- tryCatch({
     CVXR::solve(prob, solver="ECOS",
                 FEASTOL = ecos_thresh, RELTOL = ecos_thresh, ABSTOL = ecos_thresh)
  }, error=function(err){
    err$message = paste(err$message,
                        "\n", "Lasso solver using ECOS has failed." ,sep="")
    cat(err$message, fill=TRUE)
    return(NULL)
  })

  ## If anything is wrong, flag to use SCS solver
  scs = FALSE
  if(is.null(result)){
    scs = TRUE
  } else {
    if(result$status != "optimal") scs = TRUE
  }

  ## Use the SCS solver
  if(scs){
    result = CVXR::solve(prob, solver="SCS", eps = scs_eps)
    if(any(is.na(result$getValue(betamat)))){ ## A clumsy way to check.
      stop("Lasso solver using both ECOS and SCS has failed.", sep="")
    }
  }
  betahat <- result$getValue(betamat)
  return(betahat)
}



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
    y = y + rnorm(length(y), 0, 0.01)
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

