##' Trim data so that both the list of binned responses and counts don't have
##' any zeros.
##'
##' @param ybin_list List of binned data, assumed to be \code{d^3} lengthed.
##' @param counts_list Counts of binned data, assumed to be \code{d^3} lengthed.
##'
##' @return Trimmed data, each of different length.
trim <- function(ybin_list, counts_list){
  assertthat::assert_that(identical(sapply(ybin_list, nrow),
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


##' Helper to check whether ylist and countslist have been trimmed of zero,
##' e.g. using \code{trim()}.
##'
##' @inheritParams flowmix_once
##'
check_trim <- function(ylist, countslist){
  assertthat::assert_that(all.equal(sapply(ylist, nrow),
                        sapply(countslist, length))==TRUE)
  num_zeros_in_counts = sapply(countslist,function(a) sum(a==0))
  assertthat::assert_that(all(num_zeros_in_counts==0))
}


##' Prettifying the coefficients alpha and beta.
##'
##' @param alpha Alpha coefficients.
##' @param beta Beta coefficients.
##' @param numclust Number of clusters
##' @param dimdat dimension of data.
##' @param X Covariates.
##'
##' @return list containing prettified alpha and beta.
reformat_coef <- function(alpha, beta,
                          numclust, dimdat,
                          X){

  p = ncol(X)

  ## All X names
  if(is.null(colnames(X))){
    Xnames = paste0("X", 1:p)
  } else {
    Xnames = colnames(X)
    stopifnot(length(Xnames)==p)
  }

  ## Reformat betas
  beta = lapply(beta, function(b){
    colnames(b) = paste0("dim-", 1:dimdat)
    rownames(b) = c("intp", Xnames) ##paste0("X", 1:p)
    return(b)
  })
  names(beta) = paste0("clust-", 1:numclust)

  ## Reformat alphas
  if(!is.null(alpha)){
    rownames(alpha) = paste0("clust-", 1:numclust)
    colnames(alpha)[(1:(p+1))] = c("intp", Xnames)
  }

  return(list(alpha = alpha,
              beta = beta))
}


##' Inverse Value Matching; Negation of \code{%in%}. Returns the elements of
##' \code{x} that are not in \code{y}.
##'
##' @param x a vector
##' @param y a vector
##'
##' @return Elements of \code{x} that are not in \code{y}.
##' @noRd
'%ni%' <- function(x, y){
  ## return(Negate('%in%')(x))
  return( !(x %in% y) )
}


##' Helper function to lag a vector
##' @param x Numeric vector.
##' @param k Number of lags.
##'
##' @return Lag-padded numeric vector.
##'
lagpad <- function(x, k) {
  if (k>0) {
    return (c(rep(NA, k), x)[1 : length(x)] );
  }
  else {
    return (c(x[(-k+1) : length(x)], rep(NA, -k)));
  }
}

##' Symmetric difference.
##'
##' @param a One vector.
##' @param b Another vector.
##'
##' @return Symmetric difference of discrete set \code{a} and \code{b}.
##'
##' @export
##'
sym_diff <- function(a, b){
  unique(c(setdiff(a,b), setdiff(b,a)))
}


##' Check if zero pattern in coefficients, across EM iterations, has stabilized.
##'
##' @param zero.betas List of patterns in beta, over EM interations.
##' @param zero.alphas List of zero patterns in alpha, over EM iterations.
##' @param iter Iteration number.
##'
##' @return TRUE if both beta and alpha coefficients' zero patterns have
##'   stabilized.
check_zero_stabilize <- function(zero.betas, zero.alphas, iter){

  ## Check beta
  beta.sym.diffs = Map(sym_diff, zero.betas[[iter]], zero.betas[[iter-1]])
  num.beta.sym.diffs = sapply(beta.sym.diffs, length)
  zero.beta.stable = all(num.beta.sym.diffs <= 1)

  ## Check Alpha
  zero.alpha.stable = (length(sym_diff(zero.alphas[[iter]], zero.alphas[[iter-1]])) <= 1)

  return(zero.alpha.stable & zero.beta.stable)
}


##' A helper function to print the progress of a loop or simulation.
##'
##' @param isim Replicate number.
##' @param nsim Total number of replicates.
##' @param type Type of job you're running. Defaults to "simulation".
##' @param lapsetime Lapsed time, in seconds (by default).
##' @param lapsetimeunit "second".
##' @param start.time start time.
##' @param fill Whether or not to fill the line.
##'
##' @return No return
print_progress <- function(isim, nsim,
                           type = "simulation", lapsetime = NULL,
                           lapsetimeunit = "seconds", start.time = NULL,
                           fill = FALSE){

    ## If lapse time is present, then use it
    if(fill) cat(fill = TRUE)
    if(is.null(lapsetime) & is.null(start.time)){
            cat("\r", type, " ", isim, "out of", nsim)
    } else {
        if(!is.null(start.time)){
            lapsetime = round(difftime(Sys.time(), start.time,
                                       units = "secs"), 0)
            remainingtime = round(lapsetime * (nsim-isim)/isim,0)
            endtime = Sys.time() + remainingtime
        }
        cat("\r", type, " ", isim, "out of", nsim, "with lapsed time",
            lapsetime, lapsetimeunit, "and remaining time", remainingtime,
            lapsetimeunit, "and will finish at", strftime(endtime))
    }
    if(fill) cat(fill = TRUE)
}




##' Altering the beta coefficient (for all clusters) by adding zero coefficients
##' to certain rows.
##'
##' @param beta beta coefficients (list of p x dimdat matrices).
##' @param flatX Indices, out of 1 through p, of which variables are flat.
##' @param orig_names Original column names of X.
##'
##' @return Same format as beta, but all-zero rows added to beta matrix in each
##'   cluster.
alter_beta <- function(beta, flatX, orig_names){

  ## Basic checks
  stopifnot(length(flatX) >= 1)
  stopifnot(nrow(beta[[1]]) + length(flatX) == 1 + length(orig_names)) ## p+1
  dimdat = ncol(beta[[1]])
  numclust = length(beta)
  stopifnot(all(rownames(beta[[1]]) %in% c("intp", orig_names)))

  ## Make new rows
  newrows = matrix(0, ncol = dimdat, nrow = length(flatX))
  rownames(newrows) = orig_names[flatX]

  ## Add it to each beta coefficient
  for(iclust in 1:numclust){
    beta[[iclust]] = rbind(beta[[iclust]], newrows)
    beta[[iclust]] = (beta[[iclust]])[c("intp", orig_names),]
  }

  ## Return the altered beta
  return(beta)
}


##' Altering the beta coefficient (for all clusters) by adding zero coefficients
##' to certain rows.
##'
##' @param alpha alpha coefficients (One numclust x p matrix).
##' @param flatX Indices, out of 1 through p, of which variables are flat.
##' @param orig_names Original column names of X.
##'
##' @return Same format as alpha, but with zeros in the columns where variable
##'   that was flat.
alter_alpha <- function(alpha, flatX, orig_names){

  ## Basic check
  numclust = nrow(alpha)
  stopifnot(length(flatX) >= 1)
  stopifnot(all(colnames(alpha) %in% c("intp", orig_names)))
  stopifnot(ncol(alpha) + length(flatX) == length(orig_names) + 1) ## p+1

  ## Make new columns
  newcol = matrix(0, nrow = numclust, ncol = length(flatX))
  colnames(newcol) = orig_names[flatX]

  ## Add it to alpha
  alpha = cbind(alpha, newcol)
  alpha = alpha[,c("intp", orig_names)]

  ## Return the altered alpha
  return(alpha)
}
