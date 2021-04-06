context("Broad-stroke tests that the main EM algorithm is sound.")

test_that("Algorithm works on 1d, 2d and 3d data", {

  for(dimdat in 1:3){

    ## Generate 2d data if dimdat = 1
    dt = generate_data_generic(p = 5, TT = 50, fac = 1, nt = 100,
                               dimdat = dimdat + (if(dimdat==1) 1 else 0))
    ylist = dt$ylist

    ## Reduce 2d data first if dimdat = 1
    if(dimdat==1) ylist = lapply(ylist, function(y) y[,1, drop=FALSE])
    X = dt$X

    ## Fit model
    obj = flowmix_once(ylist = ylist,
                       X = X, numclust = 2,
                       prob_lambda = 0.01,
                       mean_lambda = 0.01,
                       maxdev = 0.5,
                       niter = 3,
                       verbose = FALSE)

    ## No checks, jut make sure it doesn't throw an error.
  }
})
