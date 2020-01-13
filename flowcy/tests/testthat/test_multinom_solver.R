context("Test internal sparse multinomial regression solver.")

test_that("Sparse multinomial solver helper work as expected.", {

  multinom_objective <- function(alpha, x, y, lambda,
                                          exclude.from.penalty=NULL) {
    n <- nrow(x)
    p <- ncol(x)
    L <- ncol(y)
    eta <- x %*% alpha
    v = 1:p
    if(!is.null(exclude.from.penalty)){
      stopifnot(all(exclude.from.penalty %in% (1:p)))
      v = (1:p)[-exclude.from.penalty]
   }
    ys = rowSums(y)
    (1/n) * (sum(eta * y) - sum(ys * log(rowSums(exp(eta))))) -
      lambda * sum(abs(alpha[v,]))
  }





  ## Generate some data
  TT = 20
  p = 3
  numclust = 2
  lambda = 0.1
      seed=0

  ## Generate some dummy data.
  for(seed in c(0,1,2,3,4,5,6,7,8,9,10)){
  set.seed(seed)
  print("====== Seed =============================================")
  print(seed)
  X = matrix(rnorm(p * TT), ncol=p, nrow=TT)
  X = scale(X)
  Xa = cbind(1, X)

  eta <- Xa %*% c(1, 1, 0, 0)
  prob <- 1/(1 + exp(-eta))
  yy <- runif(TT) < prob
  y <- cbind(yy, 1 - yy)
  wts = rep(c(1,200), 10)
  y2 <- y * wts

  ## CURRENT ISSUE: when y do not have row sum of 1, results differ  ##########

  ## Without intercept ########################################################

  ## Solve it two different ways
  res1 = solve_multinom(y, X, lambda, intercept=FALSE, ntlist = wts)
  res2 = cvxr_multinom(y2, X, lambda)

  ## ## These two are the same
  ## res1 = solve_multinom(y2, X, lambda, intercept=FALSE, ntlist=1/wts)
  ## res2 = cvxr_multinom(y, X, lambda)

  print("glmnet")
  print(round(res1[-1,],3))
  print("cvxr")
  print(round(res2,3))
  print("cvxr old")
  print(round(res3,3))

  ## Compare fitted coefficients
  ## expect_true(max(abs(res2 - res1[-1,])) < 1E-2)

  ## Compare objective values.
  objective1 = multinom_objective(res1[-1,], X, y, lambda)
  objective2 = multinom_objective(res2, X, y, lambda)
  objective3 = multinom_objective(res3, X, y, lambda)
  ## expect_true(abs((objective2 - objective1)/objective1) < 1E-4)
  print("glmnet")
  print(objective1)
  print("cvxr")
  print(objective2)
  print("cvxr old")
  print(objective3)

  ## With intercept ########################################################


  ## ## Temporary
  ## ntlist = sapply(ylist2, nrow)
  ## res1 = solve_multinom(X = as.matrix(X[,5:7]), y = resp.sum,
  ##                       lambda = lambda, intercept = TRUE,
  ##                       ntlist = ntlist)
  ## scaled.resp.sum = resp.sum/ntlist
  ## Xa = cbind(1, as.matrix(X[,5:7]))
  ## res2 = cvxr_multinom(resp.sum, Xa, lambda, 1)
  ## res2[-1,]
  ## res1[-1,]
  ## ## End of Temporary

  ## Solve it two different ways
  res1 = solve_multinom(y, X, lambda, intercept=TRUE)
  res2 = cvxr_multinom(y, Xa, lambda, 1)

  ## Compare fitted coefficients.
  expect_true(max(abs(res2[-1,] - res1[-1,])) < 2E-2)

  ## Compare objective values.
  objective1 = multinom_objective(res1, Xa, y, lambda, exclude.from.penalty=1)
  objective2 = multinom_objective(res2, Xa, y, lambda, exclude.from.penalty=1)
  expect_true(abs((objective2-objective1)/objective1) < 1E-3)


  ## When rowsums not 1, with intercept #####################################
  set.seed(0)
  y = matrix(runif(TT * 2), ncol=2)

  ## Solve it two different ways
  for(lambda in c(0.01) * c(2:10)){
    res1 = solve_multinom(y, X, lambda, intercept=TRUE)
    res2 = cvxr_multinom(y, Xa, lambda, 1)

    ## Compare fitted coefficients.
    expect_true(max(abs(res2[-1,] - res1[-1,])) < 1E-2)

    ## Compare objective values.
    objective1 = multinom_objective(res1, Xa, y, lambda, exclude.from.penalty=1)
    objective2 = multinom_objective(res2, Xa, y, lambda, exclude.from.penalty=1)
    expect_true(abs((objective2-objective1)/objective1) < 1E-3)
  }

  ## When rowsums not 1 and no intercept #####################################

  ## Solve it two different ways
  for(lambda in (3:10) * 0.01){
    res1 = solve_multinom(y, X, lambda, intercept=FALSE)
    res2 = cvxr_multinom(y, X, lambda)

    ## Compare fitted coefficients.
    expect_true(max(abs(res2 - res1[-1,])) < 1E-2)

    ## Compare objective values.
    objective1 = multinom_objective(res1[-1,], X, y, lambda, exclude.from.penalty=1)
    objective2 = multinom_objective(res2, X, y, lambda, exclude.from.penalty=1)
    expect_true(abs((objective2-objective1)/objective1) < 1E-3)
  }

})
