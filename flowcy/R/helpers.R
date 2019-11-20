##' Temporary helper function for printing memory
mymem <- function(msg){
  if(!is.null(msg))print(msg)
  print(gc()["Vcells","(Mb)"])
}

##' Trim data so that both the list of binned responses and counts don't have
##' any zeros.
trim <- function(ybin_list, counts_list){
  assert_that(identical(sapply(ybin_list, nrow),
                        sapply(counts_list, length)))
  for(tt in 1:TT){

    counts = counts_list[[tt]]
    ybin = ybin_list[[tt]]

    ## Trim both if necessary
    if(any(counts==0)){
      inds = which(counts==0)
      ybin_list[[tt]] = ybin[-which(counts==0), ]
      counts_list[[tt]] = as.numeric(counts[-which(counts==0)])
    }
  }
  return(list(ybin_list= ybin_list, counts_list=counts_list))
}


##' Helper to check whether ylist and countslist have been trimmed of zero
##' (e.g. using \code{trim()}).
check_trim <- function(ylist, countslist){
  assert_that(all.equal(sapply(ylist, nrow),
                        sapply(countslist, length))==TRUE)
  num_zeros_in_counts = sapply(countslist,function(a) sum(a==0))
  assert_that(all(num_zeros_in_counts==0))
}
