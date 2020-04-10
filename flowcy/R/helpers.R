##' Temporary helper function for printing memory
mymem <- function(msg){
  if(!is.null(msg))print(msg)
  print(gc()["Vcells","(Mb)"])
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


##' Prettifying the coefficients alpha and beta.
##'
##' @param alpha Alpha coefficients.
##' @param beta Beta coefficients.
##' @param p p
##' @param numclust numclust
##' @param dimdat dimdat
##'
##' @return list containing prettified alpha and beta.
reformat_coef <- function(alpha, beta,
                          p, numclust, dimdat,
                          X){
  ## p = ncol(best.res$X)
  ## numclust = best.res$numclust
  ## dimdat = best.res$dimdat


  ## All X names
  if(is.null(colnames(X))){
    Xnames = paste0("X", 1:p)
  } else {
    Xnames = colnames(X)
    stopifnot(length(Xnames)==p)
  }

  ## Reformat betas
  beta = lapply(beta, function(b){
    ## b = round(b,2)
    colnames(b) = paste0("dim-", 1:dimdat)
    rownames(b) = c("intp", Xnames) ##paste0("X", 1:p)
    return(b)
  })
  names(beta) = paste0("clust-", 1:numclust)

  ## Reformat alphas
  ## alpha = round(best.res$alpha,2)
  if(!is.null(alpha)){
    rownames(alpha) = paste0("clust-", 1:numclust)
    colnames(alpha)[(1:(p+1))] = c("intp", Xnames)
  }

  return(list(alpha = alpha,
              beta = beta))
}
