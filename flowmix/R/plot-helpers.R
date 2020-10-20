##' Collapses from a 3d cytogram to two dimensions. This is mainly used by
##' \code{scatterplot_2d()}.
##'
##' @param y 3d cytogram.
##' @param counts The multiplicity for each of the particles in \code{y}.
##' @param dims Two of \code{c(1:3)}.
##'
##' @return 3-column matrix; first two columns are the dimensions in
##'   \code{dims}.
##'
##' @import dplyr
##' @export
collapse_3d_to_2d <- function(y, counts, dims=1:2){

  ## Basic checks
  stopifnot(all(dims %in% 1:3))
  stopifnot(length(dims)==2)

  ## Aggregate
  ymat = cbind(y[,dims], counts) %>% as.data.frame()
  names(ymat)[1:2] = c("dim1", "dim2")
  ymat_summary <- ymat %>% dplyr::group_by(dim1, dim2) %>%
    dplyr::summarise(counts=sum(counts), .groups = 'drop') %>%
    as.matrix()

  ## Basic check
  if(is.null(colnames(y)) | all(colnames(y)=="")){
    colnames(ymat_summary)[1:2] = paste0("dim", dims)
  } else {
    colnames(ymat_summary)[1:2] = colnames(y)[dims]
  }

  return(ymat_summary)
}
