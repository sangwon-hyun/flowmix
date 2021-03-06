##' Bootstrapping residuals using a coin-flip assignment of cytogram particles.
##'
##' @param ylist_orig Original data, containing the actual particles.
##' @param res A \code{flowmix} class object.
##' @param countslist Counts.
##'
##' @return \code{ylist_bootstrapped} is a dataset that is exactly the same size
##'   as \code{ylist_orig}, but with the particles created by bootstrapping
##'   residuals (after having pooled all of them, across all time points).
##'
##' @examples
##' \dontrun{
##' ## Generate data and fit model
##' set.seed(0)
##' dat1 = generate_data_generic(dimdat = 2, prob1 = 1/8, nt=2000)
##' ylist = dat1$ylist
##' X = dat1$X
##' res = flowmix(ylist, X, numclust=4, niter=300,
##'               mean_lambda = 5E-3, prob_lambda = 5E-3)
##'
##' ## Now, create the bootstrapped datasets.
##' ylists <- lapply(1:nsim, function(isim){
##'   new_ylist = bootstrap(ylist, res)$ylist
##' })
##'
##' ## Make a list of bootstrapped models.
##' bootreslist <- lapply(1:nboot, function(iboot){
##'   ylist = ylists[[iboot]]
##'   res = flowmix(ylist, X, numclust=4, niter=300,
##'                 mean_lambda = 5E-3, prob_lambda = 5E-3)
##' })
##' }
##'
##' @importFrom magrittr %>%
##' @importFrom magrittr %<>%
##'
##' @export
##'
bootstrap <- function(ylist, res, countslist = NULL, verbose=FALSE){

  ## Setup
  numclust = res$numclust
  TT = length(ylist)
  dimdat = ylist %>% .[[1]] %>% ncol()

  ## Conduct the E-step once to calculate responsibilities
  resp <- Estep(res$mn, res$sigma, res$prob, ylist = ylist,
                numclust = res$numclust, first_iter = TRUE)

  ## Do a coin-flip draw of membership
  drawslist = draw_membership(resp)
  rm(resp)

  ## Draw residuals
  residuals = get_residuals(ylist, res, drawslist)

  ## Pool all the residuals
  residuals_pooled = lapply(1:res$numclust, function(iclust){
    do.call(rbind, lapply(1:TT, function(tt) residuals[[tt]][[iclust]]))
  })

  ## NEW: Obtain and pool all the biomasses (UGHH)
  if(!is.null(countslist)){

    ## Equivalent  get_residuals()
    index_by_clust <- get_index_by_clust(drawslist)

    counts_pooled = lapply(1:res$numclust, function(iclust){
      do.call(c, lapply(1:TT, function(tt){
        counts = countslist[[tt]]
        index = index_by_clust[[tt]][[iclust]]
        return(counts[index])
      }))
    })

    ## Make sure the number of pooled counts are the same as the number of
    ## pooled residuals
    assertthat::assert_that(all(sapply(counts_pooled, length) == sapply(residuals_pooled, nrow)))
  }


  ## Get the number of coin-flipped particles for each cluster
  ntklist = sapply(residuals, function(resids) resids %>% sapply(., nrow)) %>% t()

  ## Different number of clusters
  ## ntlist = ntklist %>% apply(., 1, sum)
  ## ((1/ntlist) * ntklist ) %>% matplot(type='l')
  ## (ntklist ) %>% matplot(type='l')
  ## plot(ntlist, type='l')
  ## ntlist2 = ylist %>% sapply(., nrow)

  ## Continue here!!! I think it's correct but out of memory... I should just separate out..
  new_y_and_counts_list = list()
  ## <- lapply(1:TT, function(tt){
  for(tt in 1:TT){
    if(verbose) printprogress(tt, TT, fill = TRUE)
    new_y_by_clust <- lapply(1:numclust, function(iclust){
      if(verbose) printprogress(iclust, numclust, fill = TRUE)

      ## How many particles in this cluster to pick?
      ntk = ntklist[tt, iclust]
      if(ntk == 0) return(rep(NA,dimdat) %>% rbind() %>% .[-1,])

      ## Draw the new y's and multiplicities (counts)
      bootrows = sample(1:sum(ntklist[,iclust]), ntk, replace = TRUE)
      resampled_resids <- residuals_pooled[[iclust]][bootrows,,drop=FALSE]
      new_y <- resampled_resids %>% sweep(., 2, res$mn[tt,,iclust], "+")

      ## If there are counts, include
      if(!is.null(countslist)){
        resampled_counts <- counts_pooled[[iclust]][bootrows] ## DO SOMETHING ABOUT THIS
        assertthat::assert_that(length(resampled_counts) == nrow(new_y))
        return(cbind(new_y, counts = resampled_counts))
      } else {
        return(new_y)
      }
    })
    ## return(new_y_by_clust)
    new_y_and_counts_list[[tt]] = new_y_by_clust
  }

  ## Return ylist and countslist separately
  if(!is.null(countslist)){
    new_ylist <- new_y_and_counts_list %>% lapply(., function(a){
      X = do.call(rbind, a) %>% as_tibble()
      return(X %>% select(-counts) %>% as.matrix())
    })
    new_counts <- new_y_and_counts_list %>%
      lapply(., function(a){ do.call(rbind, a)[,"counts", drop=TRUE] %>% unname() })
  } else {
    new_ylist <- new_y_and_counts_list %>% lapply(., function(a){ do.call(rbind, a)})
    new_counts = NULL
  }


  ## saveRDS(list(ylist = new_ylist,
  ##              X = res$X,
  ##              ntklist = ntklist,
  ##              drawslist = drawslist,
  ##              residuals = residuals),
  ##         file = file.path("~/Desktop/bootstrap.RDS"))

  return(list(ylist = new_ylist,
              countslist = new_counts,
              X = res$X,
              ntklist = ntklist,
              drawslist = drawslist,
              residuals = residuals))
}



##' Drawing memberships by coinflip using the responsibilities.
##'
##'
##' @param resp Posteerior cluster probabilities, or responsibilities, from
##'   \code{Estep()}.
##'
##' @importFrom magrittr %>%
##'
##' @export
draw_membership <- function(resp){

  TT = length(resp)

  drawslist = lapply(1:TT, function(tt){
    printprogress(tt, TT)
    draws = resp[[tt]] %>% apply(., 1, function(p){
     rmultinom(1, 1, p)}) %>% t()
  })

  return(drawslist)
}

##' Drawing memberships by coinflip using the "popular vote" i.e. pick
##' membership as cluster with highest responsibility..
##'
##' @importFrom magrittr %>%
##' @export
draw_membership_popular_vote <-function(resp){
  TT = length(resp)
  drawslist = lapply(1:TT, function(tt){
    draws = resp[[tt]] %>% apply(., 1, function(p){
      P = rep(0, length(p))
      P[which.max(p)] = 1
      P
    }) %>% t()
  })
  return(drawslist)
}


##' Obtain residuals of particles in \code{ylist} from the mean using model
##' \code{res}, using randomly drawn memberships stored in \code{drawslist}.
##'
##' @param ylist Data.
##' @param res Model; \code{flowmix} class object.
##' @param drawslist Membership draws, from \code{draw_membership()}.
##' @param mc.cores Number of cores, to be used by \code{parallel::mclapply()}.
##'
##' @export
get_residuals <- function(ylist, countslist = NULL, res, drawslist, mc.cores = 1){

  ## Setup
  numclust = res$numclust
  mn = res$mn
  TT = length(ylist)

  ## Obtain index from the draws
  index_by_clust <- get_index_by_clust(drawslist)

  ## Obtain all the residuals of the drawn particles.
  residuals_by_cluster = parallel::mclapply(1:TT, function(tt){
    y = ylist[[tt]]
    resid_by_clust = lapply(1:numclust, function(iclust){
      ind = index_by_clust[[tt]][[iclust]]
      this_clust_y = y[ind,, drop=FALSE]
      this_clust_residuals = sweep(this_clust_y, 2, mn[tt,,iclust])
    })
    names(resid_by_clust) = paste0("Clust", 1:numclust)
    return(resid_by_clust)
  }, mc.cores = mc.cores)

  ## Obtain the accompanying counts
  if(!is.null(countslist)){
    countslist_by_cluster = mclapply(1:TT, function(tt){
      counts = countslist[[tt]]
      counts_by_clust = lapply(1:numclust, function(iclust){
        ind = index_by_clust[[tt]][[iclust]]
        this_clust_counts = counts[ind]
      })
      names(counts_by_clust) = paste0("Clust", 1:numclust)
      return(counts_by_clust)
    }, mc.cores = mc.cores)
  } else {
    countslist_by_cluster = NULL
  }

  return(list(residuals_by_cluster = residuals_by_cluster,
              countslist_by_cluster = countslist_by_cluster))
}

##' Obtain the row numbers of each (binned) cytogram for bootstrap.
##'
##' @inheritParams get_residuals
##'
##' @export
get_index_by_clust <- function(drawslist){
  numclust = ncol(drawslist[[1]])
  TT = length(drawslist)
  index_by_cluster = lapply(1:TT, function(tt){
    draws = drawslist[[tt]]
    ind_by_clust = lapply(1:numclust, function(iclust){
      inds = which(draws[,iclust] == 1)
    })
    names(ind_by_clust) = paste0("Clust", 1:numclust)
    return(ind_by_clust)
  })
  return(index_by_cluster)
}


##' Checks whether the dataset and \code{flowmix}-class object are
##' compatible. (Under construction! Can you think of anything else?). Used in
##' unit tests.
##'
##' @param ylist Data.
##'
##' @param res Model; \code{flowmix} class object.
##'
##' @importFrom assertthat assert_that
##'
check_compatible <- function(ylist, res){
  assertthat::assert_that(res$TT == length(ylist))
  assertthat::assert_that(res$dimdat == ncol(ylist[[1]]))
}

##' Checks whether two objects of \code{flowmix} class are exactly the same
##' size.
##'
##' @param ylist1 One \code{flowmix} class object.
##'
##' @param ylist2 Another \code{flowmix} class object.
##'
##' @importFrom assertthat assert_that
##'
check_if_same_size <- function(ylist1, ylist2){

  assertthat::assert_that(length(ylist1) == length(ylist2))
  assertthat::assert_that(all(sapply(ylist1, nrow) == sapply(ylist2, nrow)))
  assertthat::assert_that(all(sapply(ylist1, ncol) == sapply(ylist2, ncol)))

}






## match_clusters <- function(newres, oldres){

##   ## Match oldres and newres
##   oldres$mn - newres$mn

##   ## Get best model
##   cvres = cv_summary(destin = file.path("~/Dropbox/research/usc/hpc-output/blockcv-2-76-5"))
##   res = cvres$bestres
##   res$prob = res$pie ## for back-compatibility

##   ## Load original dataset
##   dat = readRDS("~/repos/cruisedat/export/MGL1704-hourly-paper-1d-diam-not-binned.RDS")
##   orig_ylist = dat$ylist

##   ## Conduct the E-step once to calculate responsibilities
##   resp <- Estep(res$mn, res$sigma, res$prob, ylist = orig_ylist,
##                 numclust = res$numclust, first_iter = TRUE)
##   la('flowmix')
##   mem = draw_membership(resp)
##   mem %>% length()
##   mem %>% .[[1]] %>% dim()

##   orig_ylist[[1]]
##   myfun = function(a){which(a==1)}
##   ## for(iclust in 1:numclust){
##   TT = length(dat$ylist)
##   memlist = lapply(1:TT, function(tt){
##     onemem = lapply(1:numclust, function(iclust){
##       mem %>% .[[tt]] %>% .[,iclust] %>% myfun()
##     })
##   })
##   tt = 1
##   iclust = 1
##   orig_ylist[[tt]][Mems[[iclust]], ]

##   drawslist = draw_membership(resp)
##   draw_membership

##   ## Calculate the maximum of the problem.

## }
