##' A very rough warmstarts for covariate EM.
##' @param ylist list of data.
##' @param numclust number of clusters desired.
##' @return An array of dimension (T x dimdat x numclust).
warmstart_covar <- function(ylist, numclust){

  dimdat = ncol(ylist[[1]])
  TT = length(ylist)

  ## Collapse all the data
  all.y = do.call(rbind, ylist)

  ## Run k-means once on collapsed data.
  ## obj = kmeans(all.y, numclust)
  ## numclust = 5
  avg.num.rows = round(mean(sapply(ylist,nrow)))
  ## par(mfrow=c(5,5))
  ## for(ii in 25:1){
  ##   print(ii)
  ##   set.seed(ii)
  some.of.all.y = all.y[sample(1:nrow(all.y), avg.num.rows),]
  obj = kmeans(some.of.all.y, numclust, algorithm="MacQueen")
    ## rm(obj)


  ## ## Plot the results (temporary)
  ## plot(some.of.all.y[,1:2], type='p',cex=0.1)
  ## points(obj$centers[,1:2], col='red', pch=16)
  ##   }

  ## Repeat it TT times and return it
  centres = array(NA, dim=c(TT, dimdat, numclust))
  for(tt in 1:TT){
    centres[tt,,] = t(obj$centers)
  }
  stopifnot(dim(centres) == c(TT, dimdat, numclust)) ## Unnecessary, but still.
  return(centres)
}
