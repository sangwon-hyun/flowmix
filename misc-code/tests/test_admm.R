context("Test M-step admm solver.")

test_that("EM algorithm gives same results using ADMM solver or CVXR solver in M step.", {

  ## Generate data
  la('flowmix')
  set.seed(0)
  TT = 100
  obj = generate_data_generic(TT = TT, fac = 0.1, dimdat = 3)
  ylist = obj$ylist
  X = obj$X
  numclust = 2
  niter = 4
  for(mean_lambda in c(0,1,3,10,50,100)){
    mean_lambda=1
    print(mean_lambda)
    set.seed(1)
    res.admm = flowmix_once(ylist, X,
                            numclust = numclust,
                            niter = niter,
                            mean_lambda = mean_lambda,
                            prob_lambda = 10,
                            verbose = FALSE,
                            maxdev = 0.5,
                            admm = TRUE,
                            admm_rho = 10,
                            admm_err_rel = 1E-3)
    set.seed(1)
    res.admm.rcpp = flowmix_once(ylist, X,
                                 numclust = numclust,
                                 niter = niter,
                                 mean_lambda = mean_lambda,
                                 prob_lambda = 10,
                                 verbose = FALSE,
                                 maxdev = 0.5,
                                 admm = TRUE,
                                 admm_rho = 10,
                                 admm_err_rel = 1E-3,
                                 rcpp = TRUE)

    set.seed(1)
    res.cvxr = flowmix_once(ylist, X,
                            numclust = numclust,
                            niter = niter,
                            mean_lambda = mean_lambda,
                            prob_lambda = 10,
                            maxdev = 0.5,
                            admm = FALSE)


    for(iclust in 1:numclust){

      maxgap = max(abs(res.admm$beta[[iclust]][-1,] - res.cvxr$beta[[iclust]][-1,]))
      expect_true(maxgap < 1E-2)
    }
    maxgap = max(abs(res.admm$mn - res.cvxr$mn))
    expect_true(maxgap < 1E-1)
  }
})

## test_that("ADMM solver gives the same M-step results as the CVXR solver.", {

##   ## Main test: Compare the two beta M steps.
##   load(file.path("./test_admm.Rdata")) ## Keep this file in the test folder.
##   for(mean_lambda in c(0,1,100,200,500,1000)){
##     ## microbenchmark::microbenchmark({
##     res.beta.admm = Mstep_beta_admm(resp = resp, ylist = ylist, X = X,
##                                     mean_lambda = mean_lambda, sigma = sigma,
##                                     maxdev = 0.5,
##                                     niter = 1000,
##                                     rho = 0.1,
##                                     zerothresh = 1E-6)
##     ## }, times=10)
##     ## microbenchmark::microbenchmark({
##     res.beta.la.admm = Mstep_beta_admm(resp = resp, ylist = ylist, X = X,
##                                        mean_lambda = mean_lambda, sigma = sigma,
##                                        maxdev = 0.5,
##                                        ## rho = 100,
##                                        rho = 0.1,
##                                        zerothresh = 1E-6,
##                                        ## Inner number of iterations is niter
##                                        niter = 1000,
##                                        ## Outer number of iterations is warmstart_niter
##                                        local_adapt = TRUE,
##                                        local_adapt_niter = 20)
##     ## }, times=10)
##     res.beta.cvxr = Mstep_beta(resp = resp, ylist = ylist, X = X,
##                                mean_lambda = mean_lambda, sigma = sigma, ## maxdev = maxdev,
##                                maxdev = 0.5,
##                                sigma_eig_by_clust = sigma_eig_by_clust,
##                                zerothresh = 1E-6)

##     for(iclust in 1:numclust){
##       maxgap = max(abs(res.beta.admm$beta[[iclust]][-1,] - res.beta.cvxr$beta[[iclust]][-1,]))
##       expect_true(maxgap < 1E-2)

##       maxgap = max(abs(res.beta.la.admm$beta[[iclust]][-1,] - res.beta.cvxr$beta[[iclust]][-1,]))
##       expect_true(maxgap < 1E-2)
##     }
##     maxgap = max(abs(res.beta.admm$mns - res.beta.cvxr$mns))
##     if(mean_lambda>0) expect_true(maxgap < 3E-3) ## mean_lambda=0 gives slightly less precise results.
##     maxgap = max(abs(res.beta.la.admm$mns - res.beta.cvxr$mns))
##     if(mean_lambda>0) expect_true(maxgap < 3E-3)

##     ## ## Sparsity pattern (not checking this for now)
##     ## for(iclust in 1:numclust){
##     ##   print(iclust)
##     ##   print(table(c(which(res.beta.admm$w[[iclust]] == 0 ),
##     ##                 which(abs(res.beta.cvxr$beta[[iclust]][-1,]) < 1E-5 ))))
##     ## }

##     ## ## Plotting code (in case it's useful)
##     ## la.admm.beta <- unlist(lapply(res.beta.la.admm$beta, function(beta){beta[-1,]}))
##     ## cvxr.beta <- unlist(lapply(res.beta.cvxr$beta, function(beta){beta[-1,]}))
##     ## plot(y = cvxr.beta, x = la.admm.beta)
##     ## abline(0, 1)
##     ## plot(y = res.beta.cvxr$mns, x = res.beta.admm$mns)
##     ## abline(0, 1)
##   }
## })
