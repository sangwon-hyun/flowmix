##' A very rough warmstarts for covariate EM.
##' @param ylist list of data.
##' @param numclust number of clusters desired.
##' @return (TT x dimdat x numclust)
warmstart_covar <- function(ylist, numclust){

  dimdat = ncol(ylist[[1]])
  TT = length(ylist)

  ## Collapse all the data
  all.y = do.call(rbind, ylist)

  ## Run k-means once on collapsed data.
  obj = kmeans(all.y, numclust)

  ## Repeat it TT times and return it
  centres = array(NA, dim=c(TT, dimdat, numclust))
  for(tt in 1:TT){
    centres[tt,,] = t(obj$centers)
  }
  stopifnot(dim(centres) == c(TT, dimdat, numclust)) ## Unnecessary, but still.
  return(centres)
}
