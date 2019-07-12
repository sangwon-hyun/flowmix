##' A very rough warmstarts for covariate EM.
##' @param ylist list of data.
##' @param numclust number of clusters desired.
warmstart_covar <- function(ylist, numclust){

  dimdat = ncol(ylist[[1]])

  ## Collapse all the data
  all.y = do.call(rbind, ylist)

  ## Run k-means once on collapsed data.
  obj = kmeans(all.y, numclust)

  ## Repeat it TT times and return it
  centres = array(NA, dim=c(TT, dimdat, numclust))
  for(tt in 1:TT){
    centres[tt,,] = t(obj$centers)
  }
  return(centres)
}
