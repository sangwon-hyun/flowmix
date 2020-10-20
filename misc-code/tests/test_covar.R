context("Test EM (the main function)")

## Test for lasso and cvx lasso
test_that("flowmix() works under no noise.", {


  ## ## Test that the version with the thing has the same results
  ## la("~/repos/flowcy/flowcy")
  ## set.seed(1)
  ## dat = generate_data_generic(TT=500, p=5, nt=2000)
  ## ylist = dat$ylist
  ## X = dat$X
  ## set.seed(0)
  ## obj1 = covarem(ylist, X, numclust=4, eigenspeed=TRUE, niter=5, nrep=1, verbose=TRUE)
  ## set.seed(0)
  ## obj2 = covarem(ylist, X, numclust=4, eigenspeed=FALSE, niter=5, nrep=1, verbose=TRUE,
  ##                faster_mvn=TRUE)
  ## expect_true(max(max(obj1$alpha-obj2$alpha))<1E-4)


  ## Something like the following would go here:

  ## source("~/repos/flowcy/main/covar-artif-2d-complex.R")
  ## numclust = 3
  ## data = do.call(rbind, ylist)
  ## la("~/repos/flowcy/flowcy")
  ## ylist = lapply(ylist, function(y){
  ##   return(y[sample(1:nrow(y), 100),])
  ## })
  ## obj = generate_data_generic(TT=10, p=5)
  ## ylist = obj$ylist
  ## X = obj$X
  ## res = covarem(ylist=ylist, X=X, numclust=4, verbose=TRUE)


## pdf(file="~/Desktop/somefile.pdf")
## alldat = do.call(rbind, ylist)
## xlim = range(alldat[,1])
## ylim = range(alldat[,2])
## for(tt in 1:50){
##   plot(ylist[[tt]], ylim=ylim, xlim=xlim, pch=16, cex=0.5, col='grey50')
##   points(t(res$mn.list[[res$final.iter]][tt,,]), col='red',
##          pch=16, cex=5*res$prob.list[[res$final.iter]][tt,])
## }
## graphics.off()


})
