##' Solve the no-intercept lasso problem using CVXR. This is not a general lasso
##' function, but designed for the beta M step of the EM algorithm. The response
##' is a (dimdat T) vector, and the covariate matrix is a ((dimdat T) x (dimdat
##' p)). This solves the optimization and organizes+returns the optimization
##' variable as a ((p+1) x dimdat) matrix.
##'
##' \deqn{\min_{\beta} 1/2 \|y - X\vect(\beta)\|^2 + \lambda \|\beta\|_1}
##'
##' In addition, this includes other modifications such as the ball constraint
##' on the means over time (whose radius is handled by \code{maxdev}, and
##' \code{sel_coef} (an array the same size as the coefficients) to govern which
##' coefficients should be forced to be zero.
##'
##' Internally, it tries an ECOS solver first, then it tries a SCS solver. The
##' latter is supposedly faster (and less exact) for larger problems, but I
##' haven't seen that in our application FYI.
##'
##' @param X Covariate matrix.
##' @param y Response vector.
##' @param lambda Regularization parameter.
##' @param sel_coef Optional variable. Numeric matrix of 0's and 1's in the same
##'   size of \eqn{\beta} i.e. (p+1) x dimdat. The ones in this matrix defines
##'   the entries in the regression coefficients |beta| that are allowed to be
##'   nonzero.
##' @param maxdev Radius of ball constraint. Defaults to NULL.
##' @param refit A boolean flag for whether \code{sel_coef} should be used.
##'
##' @return Fitted beta matrix.
cvxr_lasso <- function(y, X,  lambda, Xorig=NULL,
                       exclude.from.penalty = NULL,
                       thresh = 1E-8,
                       maxdev = NULL,
                       dimdat,
                       N,
                       refit = FALSE,
                       sel_coef = NULL,
                       ecos_thresh = 1E-8,
                       scs_eps = 1E-5){

  ## Define dimensions
  n = nrow(X)
  p = ncol(X)
  TT = (n/dimdat)
  pp = p/dimdat - 1 ## The minus 1 is for the intercept

  ## Define the parameter
  betamat <- CVXR::Variable(rows=pp+1,
                            cols=dimdat)

  ## Set up exclusion from penalty, if applicable.
  if(is.null(exclude.from.penalty)){
    v = 1:p
  } else {
    assertthat::assert_that(all(exclude.from.penalty %in% (1:p)))
    v = (1:p)[-exclude.from.penalty]
  }

  ## Set up l1-penalized regression.
  obj = CVXR::sum_squares(y - X %*% CVXR::vec(betamat)) / (2 * N)
  if(!refit) obj = obj + lambda * sum(abs(CVXR::vec(betamat)[v]))

  ## Setup the Xbeta penalty.
  constraints = list()
  if(!is.null(maxdev)){
    constraints = list(CVXR::sum_entries(CVXR::square(Xorig %*% betamat[-1,]), 1) <= rep(maxdev^2, TT))
  }
  ## Forcing some coefficients to be zero, if applicable.
  if(refit){
    sel_coef = !(sel_coef) * 1
    constraints = c(constraints,
                    CVXR::vec(betamat[-1,] * sel_coef[-1,]) == rep(0, pp*dimdat))
  }

  ## Try all two CVXR solvers.
  prob <- CVXR::Problem(CVXR::Minimize(obj), constraints)
  result = NULL
  result <- tryCatch({
     solve(prob, solver="ECOS",
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
    result = solve(prob, solver="SCS", eps = scs_eps)
    if(any(is.na(result$getValue(betamat)))){ ## A clumsy way to check.
      stop("Lasso solver using both ECOS and SCS has failed.", sep="")
    }
  }
  betahat <- result$getValue(betamat)
  return(betahat)
}
