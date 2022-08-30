##' Reorder the cluster numbers for a new flowmix object \code{newres}; the best
##' permutation (reordering) is to match the original flowmix object
##' \code{origres}.
##'
##' @param newres New flowmix object to reorder.
##' @param origres Original flowmix object.
##' @param ylist_particle The particle-level data.
##' @param fac Defaults to 100, to take 1/100'th of the particles from each time point.
##' @param verbose Loud or not?
##'
##' @return Reordered res
##'
##' @export
reorder_kl <- function(newres, origres, ylist_particle, fac = 100, verbose = FALSE){

  ## Randomly sample  1/100 of the original particles (mainly for memory reasons)
  TT = length(ylist_particle)
  N = sapply(ylist_particle, nrow) %>% sum()
  ntlist = sapply(ylist_particle, nrow)
  indlist = lapply(1:TT, function(tt){
    nt = ntlist[[tt]]
    ind = sample(1:nt, round(nt / fac), replace=FALSE)
  })

  ## Sample responsibilities
  ylist_particle_small = Map(function(ind, y){ y[ind,,drop = FALSE]  }, indlist, ylist_particle)

  ## Calculate new responsibilities
  resp_orig_small <- Estep(origres$ mn, origres$sigma, origres$prob,
                           ylist = ylist_particle_small,
                           numclust = origres$numclust, first_iter = TRUE)
  resp_new_small <- Estep(newres$mn, newres$sigma, newres$prob,
                          ylist = ylist_particle_small,
                          numclust = newres$numclust, first_iter = TRUE)
  assertthat::assert_that(all(sapply(resp_orig_small, dim) == sapply(resp_new_small, dim)))

  ## Get best ordering (using symm. KL divergence and Hungarian algorithm for
  ## matching)
  best_ord <- get_best_match_from_kl(resp_new_small, resp_orig_small)

  if(verbose) cat("New order is", best_ord, fill=TRUE)
  newres_reordered_kl = newres %>% reorder_clust(ord = best_ord)

  ## Return the reordered object
  return(newres_reordered_kl)
}



##' Reorder the results of one flowmix object so that cluster 1 through
##' \code{numclust} is in a particular order. The default is decreasing order of
##' the averages (over time) of the cluster means.
##'
##' @param res flowmix object
##' @param ord Defaults to NULL. Use if you have an ordering in mind.
##'
##' @return Same object, but with clusters reordered.
##'
##' @export
reorder_clust <- function(res, ord = NULL){

  ## Find an order by sums (averages)
  if(is.null(ord)) ord = res$mn[,1,] %>% colSums() %>% order(decreasing=TRUE)
  if(!is.null(ord)) all(sort(ord) == 1:res$numclust)

  ## Reorder mean
  res$mn = res$mn[,,ord, drop=FALSE]

  ## Reorder sigma
  res$sigma = res$sigma[ord,,,drop=FALSE]

  ## Reorder prob
  res$prob = res$prob[,ord, drop=FALSE]

  ## Reorder the alpha coefficients
  res$alpha = res$alpha[ord,, drop=FALSE] ## Also rename the row names
  rownames(res$alpha) = paste0("clust-", 1:res$numclust)

  ## Reorder the beta coefficients
  res$beta = res$beta[ord]
  names(res$beta) = paste0("clust-", 1:res$numclust)

  return(res)
}

##' Compute KL divergence from responsibilities between two models'
##' responsibilities \code{resp_new} and \code{resp_old}.
##'
##' @param resp_new New responsibilities
##' @param resp_orig Original responsiblities.
##'
##' @return Calculate reordering \code{o} of the clusters in model represented
##'   by \code{resp_new}. To be clear, \code{o[i]} of new model is the best
##'   match with the i'th cluster of the original model.
##'
##' @export
##' @importFrom clue solve_LSAP
get_best_match_from_kl <- function(resp_new, resp_orig){

  ## Basic checks
  . = NULL ## Fixing check()
  assertthat::assert_that(all(sapply(resp_new, dim) == sapply(resp_orig , dim)))

  ## Row-bind all the responsibilities to make a long matrix
  distmat = form_symmetric_kl_distmat(resp_orig %>% do.call(rbind,.),
                                      resp_new %>% do.call(rbind,.))

  ## Use Hungarian algorithm to solve.
  fit <- clue::solve_LSAP(distmat)
  o <- as.numeric(fit)

  ## Return the ordering
  return(o)
}


##' From two probability matrices, form a (K x K) distance matrix of the
##' (n)-vectors. The distance between the vectors is the symmetric KL
##' divergence.
##'
##' @param mat1 Matrix 1 of size (n x K).
##' @param mat2 Matrix 2 of size (n x K).
##'
##' @return K x K matrix containing symmetric KL divergence of each column of
##'   \code{mat1} and \code{mat2}.
form_symmetric_kl_distmat <- function(mat1, mat2){

  ## Manually add some small, in case some columns are all zero
  mat1 = (mat1 + 1E-10) %>% pmin(1)
  mat2 = (mat2 + 1E-10) %>% pmin(1)

  ## Calculate and return distance matrix.
  KK1 = ncol(mat1)
  KK2 = ncol(mat2)
  distmat = matrix(NA, ncol=KK2, nrow=KK1)
  for(kk1 in 1:KK1){
    for(kk2 in 1:KK2){
      mydist = symmetric_kl(mat1[,kk1, drop=TRUE], mat2[,kk2, drop=TRUE])
      distmat[kk1, kk2] = mydist
    }
  }
  stopifnot(all(!is.na(distmat)))
  return(distmat)
}

##' Symmetric KL divergence, of two probability vectors.
##'
##' @param vec1 First probability vector.
##' @param vec2 Second prbability vector.
##'
##' @return Symmetric KL divergence (scalar).
symmetric_kl <- function(vec1, vec2){
  stopifnot(all(vec1 <= 1) & all(vec1 >= 0))
  stopifnot(all(vec2 <= 1) & all(vec2 >= 0))
  kl <- function(vec1, vec2){
    sum(vec1 * log(vec1 / vec2))
  }
  return((kl(vec1, vec2) + kl(vec2, vec1))/2)
}
