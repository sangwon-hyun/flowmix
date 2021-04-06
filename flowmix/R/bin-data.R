##' Helper to make grid. Takes 3-lengthed list of ranges (2 length vectors
##' containing min and max) and returns a 3-lengthed list of grid points.
##'
##' @param gridsize Size of equally sized grid, to be formed in each dimension
##'   of data.
##'
##' @inheritParams flowmix_once
##'
##' @export
make_grid <- function(ylist, gridsize = 5){


  ## Getting ranges over ylist
  ylist_ranges = lapply(ylist, function(mat) apply(mat, 2, range))

  ## Get overall range.
  ylist_collapsed = do.call(rbind, ylist_ranges)
  dimdat = ncol(ylist_collapsed)
  ranges = lapply(1:dimdat, function(ii) range(ylist_collapsed[,ii]))

  ## Equally space in each dimension.
  gridpoints = lapply(ranges,
                        function(x) seq(from=x[1], to=x[2], length=gridsize+1))

  return(gridpoints)
}


##' Main function for binning a cytogram (y). Only works for y whose dimension is 3.
##'
##' @param y nt by dimdat matrix.
##' @param manual.grid grid, produced using \code{make_grid()}.
##' @param qc QC value, but in the original scale, NOT in the log scale.
##'
##' @return List containing *trimmed* ybin (a x 3) and counts (a).
##'
##' @export
bin_one_cytogram <- function(y, manual.grid, qc=NULL){

  ## Basic checks
  ## assertthat::assert_that(ncol(y) == 3)
  dimdat = ncol(y)

  ## Obtain the midpoints of each box (d x d x d array)
  midpoints <- make_midpoints(manual.grid)

  ## Count the points in each of the boxes (d x d x d array)
  counts <- make_counts(y, manual.grid, qc)

  ## Aggregate all of this into a (d^3 x 4 array)
  ybin_all <- make_ybin(counts, midpoints, colnames(y), dimdat)

  ## Extract the ybin and counts
  ybin <- ybin_all[,1:dimdat, drop=FALSE]
  counts <- ybin_all[,dimdat+1]
  stopifnot(all(colnames(y) == colnames(ybin_all)[1:dimdat]))

  obj = trim_one_cytogram(ybin = ybin, counts = counts)
  sparsecounts <- Matrix::sparseVector(counts,
                                       1:length(counts),
                                       length(counts))

  return(list(ybin = obj$ybin,
              counts = obj$counts,
              sparsecounts = sparsecounts))
}


##' Main function for binning many cytograms (\code{ylist}). Only works for
##' \code{ylist} whose dimension is 3.
##'
##' @param ylist \code{TT} lengthed list of (\code{nt} by \code{dimdat})
##'   matrices.
##' @param manual.grid grid, produced using \code{make_grid()}.
##' @param qclist Biomass (Qc) of each particle. Defaults to NULL.
##'
##' @return List containing *trimmed* ybin (a x 3) and counts (a).
##'
##' @export
bin_many_cytograms <- function(ylist, manual.grid, verbose = FALSE, mc.cores = 1, qclist = NULL){

  ## Basic checks
  TT = length(ylist)
  if(verbose) cat(fill = TRUE)
  ## assertthat::assert_that(ncol(ylist[[1]]) == 3)
  dimdat = ncol(ylist[[1]])

  ## Bin each cytogram:
  reslist = parallel::mclapply(1:TT, function(tt){
  ## reslist = lapply(1:TT, function(tt){
    if(verbose & (tt %% 10 == 0 )) print_progress(tt, TT, "binning")
    bin_one_cytogram(ylist[[tt]],
                     manual.grid = manual.grid,
                     qc = qclist[[tt]])
  }, mc.cores = mc.cores)

  ## if(is.null(qclist)) browser()

  ## Gather results and return
  ybin_list = lapply(reslist, function(res) res$ybin)
  counts_list = lapply(reslist, function(res) res$counts)
  sparsecounts_list = lapply(reslist, function(res) res$sparsecounts)
  ybin_all = make_ybin(counts = NULL,  make_midpoints(manual.grid),
                       colnames(ylist[[1]]), dimdat)

  ## Name everything
  names(ybin_list) = names(ylist)
  names(counts_list) = names(ylist)
  names(sparsecounts_list) = names(ylist)

  if(verbose) cat(fill=TRUE)
  return(list(ybin_list = ybin_list,
              counts_list = counts_list,
              sparsecounts_list = sparsecounts_list,
              ybin_all = ybin_all))
}


##########################
## Internal functions ####
##########################


##' Get the midpoints in the a grid
make_midpoints <- function(grid){
  midpoints = lapply(grid, function(x){
    midpts = apply(cbind(x, lagpad(x,-1)),1,mean)
    midpts = midpts[-which(is.na(midpts))]
    return(midpts)
  })
}


make_ybin <- function(counts, midpoints, names=NULL, dimdat = 3){
  ## if(!(dimdat %in% c(2,3))){
  ##   stop("Dimension of data needs to be 2d or 3d.")
  ## }
  if(dimdat==3){
    return(make_ybin_3d(counts, midpoints, names))
  }
  if(dimdat==2){
    return(make_ybin_2d(counts, midpoints, names))
  }
  if(dimdat==1){
    return(make_ybin_1d(counts, midpoints, names))
  }
}

##' Make a (d^3 x 4) matrix of rows that look like (x,y,z,count) from the 3
##' dimensional (d x d x d) array |counts|.
##'
##' @param counts d x d x d array containing counts. If this is NULL, then dummy
##'   counts of -100 are added.
##' @param midpoints midpoints of each bin.
##' @param names Names of the data columns.
##'
##' @return (d^3 x 4) matrix of rows that look like (x,y,z,count).
make_ybin_3d <- function(counts, midpoints, names=NULL){
  gridsize = length(midpoints[[1]])
  d = gridsize
  mat = matrix(0, nrow=d^3, ncol=4)
  mm = 1
  for(ii in 1:d){
    for(jj in 1:d){
      for(kk in 1:d){
        ## Make the row c(three coordinates, count)
        if(!is.null(counts)){
          count = counts[ii,jj,kk]
        } else {
          count = -100
        }
        mat[mm,] = c(midpoints[[1]][ii],
                     midpoints[[2]][jj],
                     midpoints[[3]][kk],
                     count)
        mm = mm + 1
      }
    }
  }
  if(!is.null(names)){
    colnames(mat) = c(names, "")
  }
  return(mat)
}


##' The same as \code{make_ybin_3d()} but in 2d.
##' @inheritParams make_ybin_3d
make_ybin_2d <- function(counts, midpoints, names=NULL){
  gridsize = length(midpoints[[1]])
  d = gridsize
  dimdat = 2
  mat = matrix(0, nrow=d^dimdat, ncol=dimdat+1)
  mm = 1
  for(ii in 1:d){
    for(jj in 1:d){
        ## Make the row c(three coordinates, count)
        if(!is.null(counts)){
          count = counts[ii,jj]
        } else {
          count = -100
        }
        mat[mm,] = c(midpoints[[1]][ii],
                     midpoints[[2]][jj],
                     count)
        mm = mm + 1
    }
  }
  if(!is.null(names)){
    colnames(mat) = c(names, "")
  }
  return(mat)
}

##' The same as \code{make_ybin_3d()} but in 1d.
##' @inheritParams make_ybin_3d
make_ybin_1d <- function(counts, midpoints, names=NULL){
  gridsize = length(midpoints[[1]])
  d = gridsize
  dimdat = 1
  mat = matrix(0, nrow=d^dimdat, ncol=dimdat+1)
  mm = 1
  for(ii in 1:d){
    ## Make the row c(three coordinates, count)
    if(!is.null(counts)){
      count = counts[ii]
    } else {
      count = -100
    }
    mat[mm,] = c(midpoints[[1]][ii],
                 count)
    mm = mm + 1
  }
  if(!is.null(names)){
    colnames(mat) = c(names, "")
  }
  return(mat)
}




##' Takes a cytogram y that is a (nt x dimdat) matrix, and makes it into a 3d
##' array that contains the counts for each box.
##'
##' @param y Single cytogram.
##' @param grid Grid, created using \code{make_grid()}.
##'
##' @return All counts, as a 3-dimensional array.
make_counts <- function(y, grid, qc=NULL){

  ## Go through y_list, identify the closest box, add a count to the midpoint of
  ## the pertinent box (there is only one such box). Upweight by the carbon
  ## quantity (QC) if available.

  ## Helper: works on one row. Sees if a row is in a particular box
  identify_box <- function(grid, yrow){
    dimdat = length(yrow)
    sapply(1:dimdat, function(idim){
      max(which(unlist(yrow[idim]) >= grid[[idim]]))
    })
  }

  ## Cycle through all rows
  nt = nrow(y)
  nn = length(grid[[1]]) - 1
  dimdat = ncol(y)
  counts = array(0, dim=rep(nn,dimdat))
  for(ii in 1:nt){
    ijk = identify_box(grid, y[ii,])
    ijk = pmin(ijk, nn) ## fixing indexing for right-end edge points.
    if(is.null(qc))  to_add = 1
    if(!is.null(qc)) to_add = qc[ii]

    ##
    if(dimdat==3){
    counts[ijk[1], ijk[2], ijk[3]]= counts[ijk[1], ijk[2], ijk[3]] + to_add
    }
    if(dimdat==2){
    counts[ijk[1], ijk[2]]= counts[ijk[1], ijk[2]] + to_add
    }
    if(dimdat==1){
    counts[ijk[1]]= counts[ijk[1]] + to_add
    }
  }
  return(counts)
}



##' Trim single cytogram so that both the list of binned responses and counts
##' don't have any zeros.
##' @param ybin Binned data (assumed to be d^3 lengthed)
##' @param counts Counts of binned data (assumed to be d^3 lengthed)
trim_one_cytogram <- function(ybin, counts){
  ## Trim both if necessary
  if(any(counts==0)){
    inds = which(counts==0)
    ybin = ybin[-which(counts==0), ]
    counts = as.numeric(counts[-which(counts==0)])
  }
  return(list(ybin = ybin, counts = counts))
}
