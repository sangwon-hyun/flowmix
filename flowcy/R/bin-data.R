##' Takes a cytogram y that is a (nt x dimdat) matrix, and makes it into a 3d
##' array that contains the counts for each box.
##' @param y Single cytogram.
##' @param grid Grid, created using \code{make_grid()}.
##'
##' @return All counts, as a 3-dimensional array.
make_counts <- function(y, grid){
  ## ## Old way: ugly triple loop
  ## nn = length(grid[[1]])-1
  ## counts = array(0, dim=c(nn,nn,nn))
  ## nt = nrow(y)
  ## nnlast = 1
  ## dimdat = ncol(y)
  ## if(dimdat == 3)  k_nn_range = 1:nn
  ## if(dimdat == 2)  k_nn_range = NULL
  ## start.time = Sys.time()
  ## for(ii in 1:nn){
  ##   printprogress(ii, nn, start.time=start.time)
  ##   for(jj in 1:nn){
  ##   ## printprogress((ii-1)*nn + jj, nn^2, start.time=start.time)
  ##     for(kk in k_nn_range){
  ##       count <- sum(sapply(1:nt, function(irow){
  ##         myrow = y[irow,]
  ##         in_grid(myrow, ii, jj, kk, grid)
  ##       }))
  ##       counts[ii,jj,kk] = count
  ##     }
  ##   }
  ## }
  ## return(counts)

  ## New way: Go through y_list, identify the closest box, add a count to the
  ## midpoint of the pertinent box (there is only one).

  ## Helper: works on one row.
  identify_box <- function(grid, yrow){
    dimdat = length(yrow)
    sapply(1:dimdat, function(idim){
      max(which(yrow[idim] > grid[[idim]]))
    })
  }

  ## Cycle through all rows
  nt = nrow(y)
  nn = length(grid[[1]])-1
  counts = array(0, dim=c(nn,nn,nn))
  for(ii in 1:nt){
    ## printprogress(ii, nt)
    ijk = identify_box(grid, y[ii,])
    counts[ijk[1], ijk[2], ijk[3]]= counts[ijk[1], ijk[2], ijk[3]] + 1
  }
  return(counts)
}

##' Helper to make grid. Takes 3-lengthed list of ranges (2 length vectors
##' containing min and max) and returns a 3-lengthed list of grid points.
make_grid <- function(ranges, gridsize=5){
  gridpoints = lapply(ranges,
                     function(x) seq(from=x[1], to=x[2], length=gridsize+1))
}

##' Helper to see if the row is in the grid box indexed by (ii, jj, kk).
in_grid <- function(myrow, ii, jj, kk, grid){
  (grid[[1]][ii] <= myrow[1] ) & (myrow[1] <= grid[[1]][ii+1]) &
  (grid[[2]][jj] <= myrow[2] ) & (myrow[2] <= grid[[2]][jj+1]) &
  (grid[[3]][kk] <= myrow[3] ) & (myrow[3] <= grid[[3]][kk+1])
}


##' Get the midpoints in the a grid
make_midpoints <- function(grid){
  midpoints = lapply(grid, function(x){
    midpts = apply(cbind(x, lagpad(x,-1)),1,mean)
    midpts = midpts[-which(is.na(midpts))]
    return(midpts)
  })
}


##' Make a matrix of (x,y,z,count).
make_ybin <- function(y, gridsize, counts, midpoints){
  d = gridsize
  mat = matrix(0, nrow=d^3, ncol=4)
  mm = 1
  for(ii in 1:d){
    for(jj in 1:d){
      for(kk in 1:d){
        ## Make the c(3 coordinates, count)
        mat[mm,] = c(midpoints[[1]][ii],
                     midpoints[[2]][jj],
                     midpoints[[3]][kk],
                     counts[ii,jj,kk])
        mm = mm + 1
      }
    }
  }
  return(mat)
}
