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
##' @param
reorder_kl <- function(newres, origres, ylist_particle, fac = 100, verbose=FALSE){

  ## Randomly sample  1/100 of the original particles
  TT = length(ylist_particle)
  N = sapply(ylist_particle, nrow) %>% sum()
  ntlist = sapply(ylist_particle, nrow)
  indlist = lapply(1:TT, function(tt){
    nt = ntlist[[tt]]
    ind = sample(1:nt, round(nt / fac), replace=FALSE)
  })

  ## Sample responsibilities
  ylist_particle_small = Map(function(ind, y){ y[ind,,drop = FALSE]  }, indlist, ylist_particle)
  ## resp_new_small = Map(function(ind, resp){ resp[ind,]  }, indlist, resp_new)
  ## resp_orig_small = Map(function(ind, resp){ resp[ind,]  }, indlist, resp_orig)

  ## Calculate new responsibilities
  resp_orig_small <- Estep(origres$mn,
                           origres$sigma,
                           origres$prob,
                           ylist = ylist_particle_small,
                           numclust = origres$numclust,
                           first_iter = TRUE)
  resp_new_small <- Estep(newres$mn,
                          newres$sigma,
                          newres$prob,
                          ylist = ylist_particle_small,
                          numclust = newres$numclust,
                          first_iter = TRUE)
  assertthat::assert_that(all(sapply(resp_orig_small, dim) == sapply(resp_new_small, dim)))

  ## Reorder the clusters: using KL divergences
  matchres = kl_from_responsibilities(resp_new_small, resp_orig_small, fac = 1)
  kls = matchres$kls
  ordmat = matchres$ordmat

  ## Reorder orders using KL divergences
  best_ord = ordmat %>% .[which.min(kls),]
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
##' responsibilities \code{resp_new} and \code{resp_old}. Because this can be
##' computationally intensive, \code{fac=100} allows you to randomly choose only
##' 1/100 of the particles to use for this.
##'
##' @param resp_new new responsibilities
##' @param resp_orig original responsiblities.
##' @param fac Defaults to 100. This means only 1/100 randomly chosen particles
##'   are used for calculating responsibilities
##'
##' @return Calculate KL divergence (of a subset of points)
##' @export
kl_from_responsibilities <- function(resp_new, resp_orig, fac = 100){

  ## Basic checks
  assertthat::assert_that(all(sapply(resp_new, dim) == sapply(resp_orig , dim)))

  ## KL divergence.
  kl <- function(p1, p2){ sum(p1 * log(p1/p2)) }

  ## Get all permutations of 1:n.
  permutations <- function(n){
      if(n == 1){
          return(matrix(1))
      } else {
          sp <- permutations(n-1)
          p <- nrow(sp)
          A <- matrix(nrow=n*p,ncol=n)
          for(i in 1:n){
              A[(i-1)*p+1:p,] <- cbind(i,sp+(sp>=i))
          }
          return(A)
      }
  }

  ## Form all permutations
  numclust = resp_new %>% .[[1]] %>% ncol()
  ordmat <- permutations(numclust)

  ## Randomly sample  1/100 of the original particles
  TT = length(resp_new)
  ntlist = sapply(resp_new, nrow)
  indlist = lapply(1:TT, function(tt){
    nt = ntlist[[tt]]
    ind = sample(1:nt, round(nt/fac), replace=FALSE)
  })

  ## Sample responsibilities
  resp_new_small = Map(function(ind, resp){ resp[ind,]  }, indlist, resp_new)
  resp_orig_small = Map(function(ind, resp){ resp[ind,]  }, indlist, resp_orig)

  ## Resp
  resp_orig_small_combined <- do.call(rbind, resp_orig_small)
  resp_new_small_combined <- do.call(rbind, resp_new_small)

  ## For every permutation, compute the KL divergence between the resp matrices.
  kls = c()
  for(irow in 1:nrow(ordmat)){
    ord = ordmat[irow, , drop = TRUE]
    kls[irow] <- kl(resp_orig_small_combined,
                    resp_new_small_combined[,ord])
  }

  ## Return orders and KL divergences.
  return(list(kls = kls,
              ordmat = ordmat))
}
