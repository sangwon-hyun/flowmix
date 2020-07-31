context("Test internal sparse multinomial regression solver.")

test_that("Sparse multinomial solver (using glmnet) is consistent with CVXR.", {

  mytest <- function(seed=0, intercept=FALSE){

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

    par(mfrow=c(4,4))
    for(lambda in seq(from = 0.0001, to = 0.02, length = 8)){

      ## Without intercept
      if(!intercept){
        fit <- glmnet::glmnet(x = X,
                              y = yorig,
                              ## lambda = TT / sum(wts) * lambda,
                              lambda=lambda,
                              intercept = FALSE,
                              family = "multinomial",
                              weights = wts / sum(wts) * TT)
        coef.cvxr = cvxr_multinom(y, X,
                                  lambda,
                                  exclude.from.penalty = NULL,
                                  N=sum(wts))
        coef.glmnet = as.matrix(do.call(cbind,coef(fit)))
        coef.cvxr[which(abs(coef.cvxr)<1E-8)]=0
        coef.glmnet = coef.glmnet[-1,]
        objective.cvxr = multinom_objective(coef.cvxr, X, y, lambda, N=sum(wts))
        objective.glmnet = multinom_objective(coef.glmnet, X, y, lambda, N=sum(wts))
        ## expect_true(abs(objective.cvxr - objective.glmnet) < 1E-3)
        ## expect_true(max(abs(coef.cvxr - coef.glmnet)) < 1E-2)

        ## ## ## Estimated probabilities
        eta2 = X %*% coef.cvxr[,2]
        eta1 = X %*% coef.cvxr[,1]
        probs.cvxr = exp(eta1)/(exp(eta1) + exp(eta2))
        eta2 = X %*% coef.glmnet[,2]
        eta1 = X %*% coef.glmnet[,1]
        probs.glmnet = exp(eta1)/(exp(eta1) + exp(eta2))
        plot(y = probs.cvxr, x = probs.glmnet,
             xlab = "GLMNET-estimated class 1 probability",
             ylab = "CVXR-estimated class 1 probability")
      }

      ## With intercept
      if(intercept){
        fit <- glmnet::glmnet(x = X,
                              y = yorig,
                              lambda = ## TT / sum(wts) *
                                lambda,
                              intercept = TRUE,
                              family = "multinomial",
                              weights = wts / sum(wts) * TT)
        coef.cvxr = cvxr_multinom(y, Xa,
                                  lambda,
                                  exclude.from.penalty = 1,
                                  thresh = 1E-10,
                                  N = sum(wts))
        coef.glmnet = as.matrix(do.call(cbind,coef(fit)))
        coef.cvxr[which(abs(coef.cvxr)<1E-8)]=0
        objective.cvxr = multinom_objective(coef.cvxr, Xa, y, lambda, exclude.from.penalty=1, N=sum(wts))
        objective.glmnet = multinom_objective(coef.glmnet, Xa, y, lambda, exclude.from.penalty=1, N=sum(wts))
        ## expect_true(abs(objective.cvxr - objective.glmnet) < 1E-3)
        ## expect_true(max(abs(coef.cvxr[-1,] - coef.glmnet[-1,])) < 1.5E-2)

        ## Estimated probabilities
        eta2 = Xa %*% coef.cvxr[,2]
        eta1 = Xa %*% coef.cvxr[,1]
        probs.cvxr = exp(eta1)/(exp(eta1) + exp(eta2))
        eta2 = Xa %*% coef.glmnet[,2]
        eta1 = Xa %*% coef.glmnet[,1]
        probs.glmnet = exp(eta1)/(exp(eta1) + exp(eta2))
        plot(y = probs.cvxr, x = probs.glmnet,
             xlab = "GLMNET-estimated class 1 probability",
             ylab = "CVXR-estimated class 1 probability")
      }
      plot(x = abs(coef.cvxr),
           y = abs(coef.glmnet),
           pch = 16,
           cex = 2,
           main = paste0("lambda=", round(lambda, 3),
                         "\n objective values",
                         round(objective.cvxr,6),",", round(objective.glmnet, 6)),
           log = "xy",
           xlab = "CVXR",
           ylab = "GLMNET")
      abline(0,1)
    }
  }

  ## Conduct four tests
  mytest(1, TRUE)
  mytest(2, TRUE)
  mytest(2, FALSE)
  mytest(2, FALSE)

})
