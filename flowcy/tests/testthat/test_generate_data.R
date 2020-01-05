context("Test data generation")

## Test the parallel version of the cross validation
test_that("Data is generated as expected.", {

  ## Generate data
  TT = 100
  dimdat = 3 ## Also try dimdat = 2
  datobj = generate_data_generic(p=30, TT=TT, fac=0.1, nt=1000, dimdat=dimdat)
  ylist_collapsed = do.call(rbind, datobj$ylist)
  ranges = lapply(1:dimdat, function(idim){range(ylist_collapsed[,idim]) })

  ## Visualizing the datapoints and centers
  ## par(ask=TRUE)
  mnlist = lapply(datobj$betalist, function(beta){ datobj$Xa %*% beta })
  for(tt in 1:TT){
    y = datobj$ylist[[tt]][,1:2]
    plot(x=y[,1], y=y[,2],
         col=datobj$classlist[[tt]], pch=16, cex=1.5, xlim = ranges[[1]], ylim = ranges[[2]])
    for(iclust in 1:4){
      mns = mnlist[[iclust]][,1:2]
      points(x = mns[tt,1],
             y = mns[tt,2], col='yellow', pch=16, cex=3)
    }
  }
  ## Just visually check this, no formal tests here.
})
