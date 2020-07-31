context("Test that the refitting works properly")

## Test for lasso and cvx lasso
test_that("covarem() with refit=TRUE option gives the correct sparsity pattern.", {

  ## Generate data
  numclust = 4
  p = 5
  nt = 100
  set.seed(1)
  ## source("./covar-artif-2d-generic.R")
  TT = 10
  dat = generate_data_generic(TT=TT, nt=nt, fac=0.05)
  X = dat$X
  ylist = dat$ylist

  ## Standardize the covariates
  X = scale(X) %>% as.matrix()

  ## Fit initially.
  ## for(maxdev in list(NULL, 0.5){
  maxdev = NULL
  set.seed(1)
  res0 = covarem(ylist = ylist, X = X,
                 numclust = numclust,
                 mean_lambda = 0.01,
                 pie_lambda = 0.01,
                 maxdev = maxdev,
                 niter = 100,
                 verbose = TRUE, nrep = 3)
  sel_coef = get_sparsity_pattern(res0)

  ## Refit once more with the sparsity pattern of res0
  set.seed(1)
  res = covarem(ylist = ylist, X = X, numclust = numclust,
                mean_lambda = 0.01,
                pie_lambda = 0.01,
                refit = TRUE,
                maxdev = maxdev,
                sel_coef = sel_coef,
                niter = 100,
                verbose=TRUE, nrep = 1)

  ## Check that the sparsity pattern is equal, up to a permutation!
  for(iclust in 1:numclust){
    testthat::expect_true(sum(sapply(1:numclust, function(iclust0){
      all(get_sparsity_pattern(res0)$beta[[iclust0]] == get_sparsity_pattern(res)$beta[[iclust]])
    })) == 1)
  }
})
