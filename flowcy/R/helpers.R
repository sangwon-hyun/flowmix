##' Temporary helper function for printing memory
mymem <- function(msg){
  if(!is.null(msg))print(msg)
  print(gc()["Vcells","(Mb)"])
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

##' Trim data so that both the list of binned responses and counts don't have
##' any zeros.
##' @param ybin_list List of binned data (assumed to be d^3 lengthed)
##' @param counts_list Counts of binned data (assumed to be d^3 lengthed)
##'
##' @return Trimmed data, each of different length.
trim <- function(ybin_list, counts_list){
  assert_that(identical(sapply(ybin_list, nrow),
                        sapply(counts_list, length)))
  TT = length(ybin_list)
  for(tt in 1:TT){
    obj <- trim_one_cytogram(ybin_list[[tt]], counts_list[[tt]])
    ybin_list[[tt]] <- obj$ybin ## ybin[-which(counts==0), ]
    counts_list[[tt]] <- obj$counts ## as.numeric(counts[-which(counts==0)])
  }
  return(list(ybin_list = ybin_list,
              counts_list = counts_list))
}


##' Helper to check whether ylist and countslist have been trimmed of zero
##' (e.g. using \code{trim()}).
check_trim <- function(ylist, countslist){
  assert_that(all.equal(sapply(ylist, nrow),
                        sapply(countslist, length))==TRUE)
  num_zeros_in_counts = sapply(countslist,function(a) sum(a==0))
  assert_that(all(num_zeros_in_counts==0))
}
