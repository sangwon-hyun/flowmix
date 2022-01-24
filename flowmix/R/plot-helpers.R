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
collapse_3d_to_2d <- function(y, counts, dims = 1:2){

  ## Basic checks
  stopifnot(all(dims %in% 1:3))
  stopifnot(length(dims) == 2)

  ## Aggregate
  ymat = cbind(y[,dims], counts) %>% as.data.frame()
  names(ymat)[1:2] = c("dim1", "dim2")
  ymat_summary <- ymat %>% dplyr::group_by(dim1, dim2) %>%
    dplyr::summarise(counts=sum(counts), .groups = 'drop') %>%
    as.matrix()
  ## TODO fix check() messages about "no visible global variable for dim1
  ## dim2" using
  ## https://community.rstudio.com/t/data-pronoun-and-no-visible-binding-for-global-variable/85362

  ## Basic check
  if(is.null(colnames(y)) | all(colnames(y)=="")){
    colnames(ymat_summary)[1:2] = paste0("dim", dims)
  } else {
    colnames(ymat_summary)[1:2] = colnames(y)[dims]
  }

  return(ymat_summary)
}


##' Collapse 3d data to 1d data.
##'
##' @param ylist ylist
##' @param countslist countslist
##' @param idim Which dimension to collapse to (out of the
##'   \code{ncol(ylist[[1]])} dimensions.) Alternatively, this can be a column
##'   name of the tables in \code{ylist}
##'
##' @return Collapsed list of ylist and countslist
##'
##' @export
collapse_3d_to_1d <- function(ylist, countslist, idim){

  ## Basic checks
  stopifnot(length(idim) == 1)
  if(is.numeric(idim)){
    if(idim == as.integer(idim)){
      assertthat::assert_that(idim > 0 & idim <= ncol(ylist[[1]]))
    }
  } else if (assertthat::is.string(idim)){
    assertthat::assert_that(idim %in% colnames(ylist[[1]]))
  } else {
    stop("Check the |idim| parameter!")
  }

  ## Collapse all dimensions
  TT = length(ylist)
  countslist_1d = parallel::mclapply(1:TT, function(tt){
    print_progress(tt, TT)

    ## Collapse the data (only one time point, for now)
    y = ylist[[tt]]
    cc = countslist[[tt]]

    ## Collapse the counts by every unique value of y[,idim]
    lv = as.numeric(levels(factor(y[,idim, drop=TRUE])))
    stopifnot(all(sort(lv) == lv)) ## making sure that unique diam values are in
                                   ## order.
    ccnew = sapply(1:length(lv), function(ii){
      sum(cc[which(as.numeric(factor(y[,idim, drop=TRUE]))==ii)])
    })
    names(ccnew) = lv
    return(ccnew)
  }, mc.cores = 1)

  ## The coordinates are the names of the elements in |countslist|
  ylist_1d = lapply(countslist_1d, function(a) as.numeric(names(a)))

  ## Final checks
  assertthat::assert_that(all(sapply(ylist_1d, length) ==
                              sapply(countslist_1d, length)))

  ## Return the ylist
  return(list(ylist = ylist_1d,
              countslist = countslist_1d))
}



##' Helper: to get ranges for each dimension, of the list of cytograms |ylist|.
##'
##' @param ylist List of cytogram data.
##'
##' @return Matrix whose rows contain the range (min & max) of each dimension.
get_range_from_ylist <- function(ylist){
  maxs = do.call(rbind, lapply(ylist, function(a) apply(a,2,max)))
  mins = do.call(rbind, lapply(ylist, function(a) apply(a,2,min)))
  maxs = apply(maxs, 2, max)
  mins = apply(mins, 2, min)
  return(rbind(mins, maxs))
}
