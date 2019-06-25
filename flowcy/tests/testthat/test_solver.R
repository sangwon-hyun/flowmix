context("Test internal sparse regression solvers.")

## Test for lasso and cvx lasso
test_that("Lasso solver helpers work as expected.", {

  ## Lasso objective function
  lasso_objective <- function(beta, x, y, lambda, exclude.from.penalty=NULL) {

    n <- nrow(x)
    p <- ncol(x)
    if(is.null(exclude.from.penalty)){
      v = 1:p
    } else {
      assert_that(all(exclude.from.penalty %in% (1:p)))
      v = (1:p)[-exclude.from.penalty]
    }
    sumsq = function(a){sum(a*a)}
    obj <- sumsq(y - x %*% beta) / (2 * n) + lambda * sum(abs(beta[v]))
    return(obj)
  }



  ## Check match between our internal solver (glmnet_lasso) and cvx_lasso,
  ## which are solving the no-intercept problem free of any scale problems.

  n <- 100
  p <- 20
  for(seed in c(1,2,3)){
    for(lambda in c(0.01, 0.02, 0.03)){

      seed=1
      lambda=0.01

      ## Generate data
      set.seed(seed)
      x <- matrix(rnorm(n*p),n,p)
      y <- runif(n)
      xa = cbind(1,x)

      ## No-intercept version is in agreement with CVX.
      res1 = cvxr_lasso(y, x, lambda)
      res2 = glmnet::glmnet(y=y, x=x, lambda=lambda, intercept=FALSE, standardize=FALSE)
      res3 = solve_lasso(y=y, x=x, lambda=lambda, intercept=FALSE)
      mat = cbind(cvx=res1, glmnet=as.numeric(coef(res2)[-1]),
                  wrapper=c(res3$b0, res3$b))
      ## print(round(mat,4))
      testthat::expect_true(max(abs(mat[,1] - mat[,2])) < 1E-4)
      testthat::expect_true(max(abs(mat[,2] - mat[,3])) < 1E-4)

      ## Intercept version is in overall agreement with CVX.
      res1 = cvxr_lasso(y, xa, lambda, exclude.from.penalty=1)
      res2 = glmnet::glmnet(y=y, x=x, lambda=lambda, intercept=TRUE, standardize=FALSE)
      res3 = solve_lasso(y=y, x=x, lambda=lambda, intercept=TRUE)
      mat = cbind(cvxr=res1, glmnet=as.numeric(coef(res2)), wrapper=c(res3$b0, res3$b))
      ## print(round(mat,4))
      testthat::expect_true(max(abs(mat[,1] - mat[,2])) < 1E-4)
      testthat::expect_true(max(abs(mat[,2] - mat[,3])) < 1E-4)

      ## What if we asked for intercept?
      ex = c(1,3,5)
      res1 = cvxr_lasso(y, x, lambda, ex)
      res2 = solve_lasso(y=y, x=x, lambda=lambda, intercept=FALSE,
                         exclude.from.penalty=ex)
      mat = cbind(cvx=as.numeric(res1),
                  wrapper=c(res2$b))
      ## print(round(mat,4))
      testthat::expect_true(max(abs(mat[,1] - mat[,2])) < 1E-2)
      obj1 = lasso_objective(mat[,1], x, y, lambda, exclude.from.penalty=ex)
      obj2 = lasso_objective(mat[,2], x, y, lambda, exclude.from.penalty=ex)
      stopifnot(abs(obj1 - obj2) < 1E-4)
    }
  }

  ## Generate data once more
  set.seed(1234)
  n <- 100
  p <- 20
  x <- matrix(rnorm(n*p),n,p)
  y <- runif(n)

  ## Variable to exclude
  ex = c(1,3,5,8,12,14)
  pen = rep(1, p)
  pen[ex] = 0
  pen = pen / sum(pen) * (p)

  ## Test across a range of lambdas (res3 is using a wrong, unadjusted lambda)
  lambdas=seq(from=0.01,to=0.1, length=100)
  matlist = objlist = list()
  for(ii in 1:100){
    lambda = lambdas[ii]
    res1 = cvxr_lasso(y, x, lambda, ex)
    lam = lambda / pen[2]
    res2 = solve_lasso(y, x, lambda, FALSE, ex)
    ## res2 = glmnet::glmnet(y=y, x=x, lambda=lam, intercept=FALSE,
    ##               penalty.factor=pen, thresh=1E-30)
    ## res3 = glmnet::glmnet(y=y, x=x, lambda=lambda, intercept=FALSE, penalty.factor=pen, thresh=1E-30)

    mat = cbind(cvx=as.numeric(res1),
                wrapper=res2$b)
              ##   glmnet=as.numeric(coef(res2)[-1]),
              ## glmnet2=as.numeric(coef(res3)[-1]))

    obj1 = lasso_objective(mat[,1], x, y, lambda, exclude.from.penalty=ex)
    obj2 = lasso_objective(mat[,2], x, y, lambda, exclude.from.penalty=ex)
    ## obj3 = lasso_objective(mat[,3], x, y, lambda, exclude.from.penalty=ex)

    matlist[[ii]] = mat
    objlist[[ii]] = c(obj1, obj2)##, obj3)
  }

  ## Objectives should match very well
  objmat = do.call(rbind, objlist)
  expect_true(max(abs(objmat[,2]-objmat[,1]))<1E-4)

  ## Fitted values should also match very well
  errs12 = errs13 = c()
  for(ii in 1:100){
    mat = matlist[[ii]]
    testthat::expect_true(max(abs(mat[,1] - mat[,2])) < 1E-2)
    errs12[ii] = max(abs(mat[,1] - mat[,2]))
    ## errs13[ii] = max(abs(mat[,1] - mat[,3]))
  }
  expect_true(max(errs12) < 1E-2) ## errs13 should be a good deal higher than
                                  ## errs12, giving weight to using sqrt(lam)
})


test_that("Sparse multinomial solver helper work as expected.", {

  ## Calculates the objective value
  multinom_objective <- function(alpha, x, y, lambda, exclude.from.penalty=NULL) {

    n <- nrow(x)
    p <- ncol(x)
    L <- ncol(y)
    eta <- x %*% alpha

    v = 1:p
    if(!is.null(exclude.from.penalty)){
      stopifnot(all(exclude.from.penalty %in% (1:p)))
      v = (1:p)[-exclude.from.penalty]
    }

    sum(eta * y) / n - sum(log(rowSums(exp(eta)))) / n - lambda * sum(abs(alpha[v,]))
  }

  ## Generate some data
  TT = 20
  p = 3
  numclust = 2
  lambda = 0.1

  ## Generate some dummy data.
  set.seed(0)
  X = matrix(rnorm(p * TT), ncol=p, nrow=TT)
  Xa = cbind(1, X)

  ## eta <- X %*% c(1, 0, 0)
  eta <- Xa %*% c(1, 1, 0, 0)
  prob <- 1/(1 + exp(-eta))
  yy <- runif(TT) < prob
  y <- cbind(yy, 1 - yy)

  ## Without intercept:

  ## Solve it two different ways
  res1 = solve_multinom(y, X, lambda, intercept=FALSE)
  res2 = cvxr_multinom(y, X, lambda)

  ## Compare fitted coefficients
  expect_true(max(res2 - res1[-1,]) < 1E-2)

  ## Compare objective values.
  expect_true(max(multinom_objective(res1[-1,], X, y, lambda) -
                  multinom_objective(res2, X, y, lambda) ) < 1E-5)

  ## With intercept:

  ## Solve it two different ways
  res1 = solve_multinom(y, X, lambda, intercept=TRUE)
  res2 = cvxr_multinom(y, Xa, lambda, 1)

  ## Compare fitted coefficients.
  expect_true(max(res2 - res1) < 1E-2)

  ## Compare objective values.
  expect_true(max(multinom_objective(res1[-1,], X, y, lambda) -
                  multinom_objective(res2[-1,], X, y, lambda)) < 1E-5)
})
