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
  beta <- CVXR::Variable(p)
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
  obj <- CVXR::sum_squares(y - X %*% beta) / (2 * n) + lambda * sum(abs(beta[v]))
  prob <- CVXR::Problem(CVXR::Minimize(obj), constraints)
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




##' (New version) Solving the no-intercept Lasso problem using CVXR.
##' \deqn{\min_{\beta} 1/2n \|y - X\beta\|^2 + \lambda \|\beta\|_1}
##' @param X Covariate matrix.
##' @param y Response vector.
##' @param lambda regularization problem.
##' @param sel_coef The ones in this matrix defines the entries in the
##'   regression coefficients |beta| that are allowed to be nonzero.
##' @return Fitted beta.
cvxr_lasso_newer <- function(y, X, Xorig, lambda,
                             exclude.from.penalty=NULL,
                             thresh=1E-8,
                             maxdev=NULL,
                             dimdat=NULL,
                             numclust=NULL,
                             refit=FALSE,
                             sel_coef=NULL
                           ){

  ## Define dimensions
  n = nrow(X)
  p = ncol(X)
  TT = (n/dimdat)
  pp = p/dimdat - 1 ## The minus 1 is for the intercept

  ## Define the parameter
  betamat <- CVXR::Variable(rows=pp+1,
                   cols=dimdat)

  ## Define the squared loss (the main part)
  loss <- sum((y - X %*% CVXR::vec(betamat))^2) / (2 * n)

  ## Set up exclusion from penalty, if applicable.
  if(is.null(exclude.from.penalty)){
    v = 1:p
  } else {
    assert_that(all(exclude.from.penalty %in% (1:p)))
    v = (1:p)[-exclude.from.penalty]
  }

  ## Perform elastic-net regression.
  obj = CVXR::sum_squares(y - X %*% CVXR::vec(betamat)) / (2 * n)
  if(!refit) obj = obj + lambda * sum(abs(CVXR::vec(betamat)[v]))

  ## Setup the Xbeta penalty.
  constraints = list()
  if(!is.null(maxdev)){
    ## mymat = cbind(rep(1,dimdat))
    ## constraints = list(square(Xorig %*% betamat[-1,]) %*% mymat <= rep(maxdev^2,TT))
    constraints = list(CVXR::sum_entries(CVXR::square(Xorig %*% betamat[-1,]), 1) <= rep(maxdev^2,TT))
  }

  if(refit){
    sel_coef = !(sel_coef) * 1
    constraints = c(constraints,
                    CVXR::vec(betamat[-1,] * sel_coef[-1,]) == rep(0, pp*dimdat))
  }

  ## Solve the problem.
  prob <- CVXR::Problem(CVXR::Minimize(obj), constraints)
  ## browser()

  ## Tempoarary
  ## ## Cehck objective values.
  ## eps.list = c(1E-10, 1E-8, 1E-7, 1E-6, 1E-5, seq(from=1E-3, to=1E-2, length=5))
  ## reslist = lapply(eps.list, function(eps){
  ##   print(eps)
  ##   start.time = Sys.time()
  ##   set.seed(0)
  ##   result <- solve(prob, solver = "SCS",
  ##                    eps = eps)
  ##   solution <- result$getValue(betamat)
  ##   objective_value = result$value
  ##   lapse.time = round(difftime(Sys.time(), start.time,
  ##                              units = "secs"), 2)
  ##   return(list(lapse.time = lapse.time,
  ##               objective_value = objective_value,
  ##               solution = solution))
  ## })

  ## eps.list = c(1E-10, 1E-8, 1E-7, 1E-6, 1E-5, seq(from=1E-3, to=1E-2, length=5))
  ## reslist_ecos = lapply(eps.list, function(eps){
  ##   print(eps)
  ##   start.time = Sys.time()
  ##   set.seed(0)
  ##   result <- solve(prob, solver="ECOS",
  ##                   FEASTOL = eps, RELTOL = eps, ABSTOL = eps)
  ##   solution <- result$getValue(betamat)
  ##   objective_value = result$value
  ##   lapse.time = round(difftime(Sys.time(), start.time,
  ##                              units = "secs"), 2)
  ##   return(list(lapse.time = lapse.time,
  ##               objective_value = objective_value,
  ##               solution = solution))
  ## })

  ## pdf("~/Desktop/ecos-by-tolerance.pdf", width=12, height=4)
  ## par(mfrow=c(1,3))

  ## lapsetimes = sapply(reslist_ecos, function(a) a$lapse.time)
  ## plot(y = lapsetimes, x = eps.list, type='o', ylim = c(0, max(lapsetimes)), main="Solve times (sec)", lwd=2, log="x", xlab = "Convergence tolerance for ECOS solver")
  ## abline(h=0, lwd=2, col='grey')
  ## abline(v=1E-8, log="x", col='pink')

  ## objectives = sapply(reslist_ecos, function(a) a$objective_value)
  ## plot(y = objectives, x = eps.list, type='o',
  ##      main="Objective value (minimization)", lwd=2, log="x",xlab = "Convergence tolerance for ECOS solver")
  ## abline(v=1E-8, log="x", col='pink')

  ## matplot(y = sapply(2:20, function(ii)(sapply(reslist_ecos, function(a) a$solution[ii, 2]))),
  ##         x = eps.list,
  ##         type='o', pch=16, lty=1, lwd=2, log="x", xlab = "Convergence tolerance for ECOS solver",
  ##         main = "Some fitted coefficient values",
  ##         ylab = "Values")
  ## abline(v=1E-8, log="x", col='pink')
  ## graphics.off()


  ## ## reslist = reslist2
  ## ## eps.list0 = c(eps.list[1:4], eps.list2[-1])
  ## ## reslist0 = c(reslist[1:4], reslist2[-1])
  ## pdf("~/Desktop/scs-by-tolerance.pdf", width=12, height=4)
  ## par(mfrow=c(1,3))

  ## lapsetimes = sapply(reslist, function(a) a$lapse.time)
  ## plot(y = lapsetimes, x = eps.list, type='o', ylim = c(0, max(lapsetimes)), main="Solve times (sec)", lwd=2, log="x", xlab = "Convergence tolerance for SCS solver")
  ## abline(h=0, lwd=2, col='grey')
  ## abline(v=1E-5, log="x", col='pink')

  ## objectives = sapply(reslist, function(a) a$objective_value)
  ## plot(y = objectives, x = eps.list, type='o',
  ##      main="Objective value (minimization)", lwd=2, log="x",xlab = "Convergence tolerance for SCS solver")
  ## abline(v=1E-5, log="x", col='pink')

  ## matplot(y = sapply(2:20, function(ii)(sapply(reslist, function(a) a$solution[ii, 2]))),
  ##         x = eps.list,
  ##         type='o', pch=16, lty=1, lwd=2, log="x", xlab = "Convergence tolerance for SCS solver",
  ##         main = "Some fitted coefficient values",
  ##         ylab = "Values")
  ## abline(v=1E-5, log="x", col='pink')
  ## graphics.off()

  ## End of tempoarary
  result <- solve(prob, solver="ECOS",
                  FEASTOL = thresh, RELTOL = thresh, ABSTOL = thresh)
  betahat <- result$getValue(betamat)

  return(betahat)
}
