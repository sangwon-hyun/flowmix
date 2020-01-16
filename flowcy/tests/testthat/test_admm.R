context("Test M-step admm solver.")

test_that("EM algorithm gives same results using ADMM solver or CVXR solver in M step.", {

  ## Generate data
  set.seed(0)
  TT = 100
  obj = generate_data_generic(TT = TT, fac = 0.1, dimdat = 3)
  ylist = obj$ylist
  X = obj$X

  ## ## Commented out plotting code, for now.
  ## pdf("~/Desktop/somefigure.pdf", width=30, height=30)
  ## numclust = 10
  ## par(mfrow = c(6, numclust))
  ## Make all comparisons
  niter = 4
  for(mean_lambda in c(0,1,3,10,50,100)){
    set.seed(1)
    res.admm = covarem_once(ylist, X, numclust = numclust, niter = niter,
                            mean_lambda = mean_lambda,
                            pie_lambda = 10,
                            verbose = TRUE,
                            maxdev = 0.5,
                            admm=TRUE,
                            admm_rho = 10,
                            admm_err_rel = 1E-3)

    set.seed(1)
    res.cvxr = covarem_once(ylist, X, numclust = numclust, niter = niter,
                            mean_lambda = mean_lambda,
                            pie_lambda = 10,
                            verbose = TRUE,
                            maxdev = 0.5,
                            admm = FALSE)


    for(iclust in 1:numclust){
      maxgap = max(abs(res.admm$beta[[iclust]][-1,] - res.cvxr$beta[[iclust]][-1,]))

      print('max gap')
      print(maxgap)
      expect_true(maxgap < 1E-2)
    }
    maxgap = max(abs(res.admm$mn - res.cvxr$mn))
    print('mean max gap')
    print(maxgap)
    expect_true(maxgap < 1E-1)

    ## ## Commented out plotting code, for now:
    ## for(iclust in 1:numclust){
    ## plot(x = res.cvxr$mn,
    ##      y = res.admm$mn)
    ##   abline(0,1)
    ##   plot(x = res.cvxr$beta[[iclust]],
    ##        y = res.admm$beta[[iclust]])
    ##   abline(0,1)
    ## }
  }
  ## graphics.off()
})

test_that("ADMM solver gives the same M-step results as the CVXR solver.", {

  ## Main test: Compare the two beta M steps.
  load(file.path("./test_admm.Rdata")) ## Keep this file in the test folder.
  for(mean_lambda in c(0,1,100,200,500,1000)){
    microbenchmark::microbenchmark({
    res.beta.admm = Mstep_beta_admm(resp=resp, ylist=ylist, X=X,
                                    mean_lambda = mean_lambda, sigma=sigma,
                                    maxdev = 0.5,
                                    niter = 1000,
                                    rho = 100,
                                    zerothresh = 1E-6)
    }, times=10)
    microbenchmark::microbenchmark({
    res.beta.la.admm = Mstep_beta_admm(resp=resp, ylist=ylist, X=X,
                                       mean_lambda = mean_lambda, sigma=sigma,
                                       maxdev = 0.5,
                                       rho = 100,
                                       zerothresh = 1E-6,
                                       ## Inner number of iterations is niter
                                       niter = 100,
                                       ## Outer number of iterations is warmstart_niter
                                       local_adapt = TRUE,
                                       local_adapt_niter = 5)
    }, times=10)
    res.beta.cvxr = Mstep_beta(resp=resp, ylist=ylist, X=X,
                               mean_lambda = mean_lambda, sigma = sigma, ## maxdev = maxdev,
                               maxdev = 0.5,
                               sigma_eig_by_clust = sigma_eig_by_clust,
                               zerothresh = 1E-6)

    for(iclust in 1:numclust){
      maxgap = max(abs(res.beta.admm$beta[[iclust]][-1,] - res.beta.cvxr$beta[[iclust]][-1,]))
      expect_true(maxgap < 1E-2)

      maxgap = max(abs(res.beta.la.admm$beta[[iclust]][-1,] - res.beta.cvxr$beta[[iclust]][-1,]))
      expect_true(maxgap < 1E-2)
    }
    maxgap = max(abs(res.beta.admm$mns - res.beta.cvxr$mns))
    expect_true(maxgap < 3E-3)
    maxgap = max(abs(res.beta.la.admm$mns - res.beta.cvxr$mns))
    expect_true(maxgap < 2E-3)

    ## ## Sparsity pattern (not checking this for now)
    ## for(iclust in 1:numclust){
    ##   print(iclust)
    ##   print(table(c(which(res.beta.admm$w[[iclust]] == 0 ),
    ##                 which(abs(res.beta.cvxr$beta[[iclust]][-1,]) < 1E-5 ))))
    ## }

    ## ## Plotting code (in case it's useful)
    ## la.admm.beta <- unlist(lapply(res.beta.la.admm$beta, function(beta){beta[-1,]}))
    ## cvxr.beta <- unlist(lapply(res.beta.cvxr$beta, function(beta){beta[-1,]}))
    ## plot(y = cvxr.beta, x = la.admm.beta)
    ## abline(0, 1)
    ## plot(y = res.beta.cvxr$mns, x = res.beta.admm$mns)
    ## abline(0, 1)
  }
})
