##' Drawing memberships by coinflip using the responsibilities.
##'
##'
##' @param resp Posteerior cluster probabilities, or responsibilities, from
##'   \code{Estep()}.
##' @param verbose Loud or not?
##'
##' @return List containing cluster membership draws.
##'
##' @export
draw_membership <- function(resp, verbose = FALSE){

  . = NULL ## Fixing check()
  TT = length(resp)

  drawslist = lapply(1:TT, function(tt) c())
  for(tt in 1:TT){
    if(verbose) print_progress(tt, TT, "Membership draw, time point = ")
    if(nrow(resp[[tt]]) == 0){
      next
    } else {
      draws = resp[[tt]] %>% apply(., 1, function(p){stats::rmultinom(1, 1, p)}) %>% t()
      drawslist[[tt]] = draws
    }
  }
  if(verbose)cat(fill = TRUE)

  return(drawslist)
}



##' Obtain residuals of particles in \code{ylist} from the mean using model
##' \code{res}, using randomly drawn memberships stored in \code{drawslist}.
##'
##' @param ylist Data.
##' @param res Model; \code{flowmix} class object.
##' @param drawslist Membership draws, from \code{draw_membership()}.
##' @param mc.cores Number of cores, to be used by \code{parallel::mclapply()}.
##' @param countslist Multiplicity for particles in \code{ylist}.
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
    countslist_by_cluster = parallel::mclapply(1:TT, function(tt){
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
##' @inheritParams get_residuals
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





##' From subsampling bootstrap results (summary files), produce stability
##' measures.
##'
##' @param nsim Number of simulations.
##' @param outputdir Contains files named "summary-(isim).RDS".
##' @param origres Model estimated from the entire dataset.
##' @param ylist_particle Original particle-level cytogram data.
##' @param X Accompanying covariate data.
##' @param prefix The prefix of the summary files, e.g. "summary-sim-" for file that are named "summary-sim-95.RDS".
##'
##' @return List of nonzero frequencies, and alpha and beta coefficients of
##'   reordered models.
##'
##' @export
get_frequency <- function(nsim, outputdir, origres, ylist_particle, X, prefix="summary-sim-"){

  ## Setup
  . = NULL ## Fixing check()
  numclust = origres$numclust

  ######################
  ###  Frequencies #####
  ######################
  beta_list = alpha_list = mn_list = list()
  start.time = Sys.time()
  for(isim in 1:nsim){
    print_progress(isim, nsim, start.time=start.time)

    ## Load the data.
    ## resfile = file.path(outputdir, paste0("summary-", isim, ".RDS"))
    filename = paste0(prefix, isim, ".RDS")
    resfile = file.path(outputdir, filename)
    if(!file.exists(resfile)) next
    cvres = readRDS(file = resfile)

    ## Estimated model
    newres = predict.flowmix(cvres$bestres, newx = X)
    ## class(newres) = "flowmix"

    ## Reorder the new res.
    newres = newres %>% reorder_kl(origres, ylist_particle, fac = 100, verbose = FALSE)

    ## Store the beta and alpha coefficients
    alpha_list[[isim]] = newres$alpha
    beta_list[[isim]] = newres$beta
    mn_list[[isim]] = newres$mn
  }

  #####################################################################
  ## Obtain the frequencies (same-sized objects as alpha and beta) ####
  #####################################################################
  threshold <- function(x, tol){abs(x) > tol}
  alpha_nonzero = lapply(alpha_list, threshold, 1E-8) %>% Reduce("+", .)

  ## Combine all betas into (K x (p+1) x d) arrays
  beta_arrays = lapply(beta_list, abind, along = 0)

  ## Calculate frequencies
  beta_freq = beta_arrays %>% purrr::map(. %>% threshold(1E-8)) %>% Reduce("+", .) %>% `/`(nsim)

  ## Combine all alphas
  alpha_freq = alpha_list %>% purrr::map(. %>% threshold(1E-8)) %>% Reduce("+", .) %>% t() %>% `/`(nsim)

  return(list(alpha_freq = alpha_freq,
              beta_freq = beta_freq,
              alpha_list = alpha_list,
              beta_list = beta_list,
              mn_list = mn_list))
}
