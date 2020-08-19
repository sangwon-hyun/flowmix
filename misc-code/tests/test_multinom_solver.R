context("Test internal sparse multinomial regression solver.")

test_that("Sparse multinomial solver (using glmnet) is consistent with CVXR.", {

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


  mytest <- function(seed){

    ## Generate some data
    TT = 1000
    p = 10
    set.seed(seed)
    X = matrix(rnorm(p * TT), ncol=p, nrow=TT)
    X = scale(X)
    Xa = cbind(1, X)

    eta <- Xa %*% c(2, 1, rep(0, p-1))
    prob <- 1/(1 + exp(-eta))
    yy <- ( runif(TT) < prob )
    yorig <- cbind(yy, 1 - yy)
    wts = rep(c(1, 10), TT / 2)
    y <- yorig * wts
    N = sum(y)

    par(mfrow=c(4,4))
    lambdas = seq(from = 0.0001, to = 0.05, length = 8)
    for(lambda in lambdas){

        ## What we use.
        coef.glmnet <- solve_multinom(y, X, lambda)
        ## What we test against.
        coef.cvxr = cvxr_multinom(y, Xa,
                                  lambda,
                                  exclude.from.penalty = 1,
                                  thresh = 1E-10,
                                  N = N)
        objective.cvxr = multinom_objective(coef.cvxr, Xa, y, lambda, exclude.from.penalty=1)
        objective.glmnet = multinom_objective(coef.glmnet, Xa, y, lambda, exclude.from.penalty=1)

        ## Estimated probabilities
        eta2 = Xa %*% coef.cvxr[,2]
        eta1 = Xa %*% coef.cvxr[,1]
        probs.cvxr = exp(eta1)/(exp(eta1) + exp(eta2))

        eta1 = Xa %*% coef.glmnet[,1]
        eta2 = Xa %*% coef.glmnet[,2]
        probs.glmnet = exp(eta1)/(exp(eta1) + exp(eta2))


      ## Make plots of estimated probabilities, and coefficients
        plot(y = probs.cvxr, x = probs.glmnet,
             xlab = "GLMNET-estimated class 1 probability",
             ylab = "CVXR-estimated class 1 probability")
      plot(x = abs(coef.cvxr),
           y = abs(coef.glmnet),
           pch = 16,
           cex = 2,
           main = paste0("lambda=", round(lambda, 3),
                         "\n objective values",
                         round(objective.cvxr,3),",", round(objective.glmnet, 3)),
           log = "xy",
           xlab = "CVXR",
           ylab = "GLMNET")
      abline(0,1)
    }
  }

  ## Condut the tests
  mytest(1)
  mytest(2)

})
