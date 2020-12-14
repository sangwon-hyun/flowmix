##' Bootstrapping residuals using a coin-flip assignment of cytogram particles.
##'
##' @param ylist_orig Original data, containing the actual particles
##' @param res A \code{flowmix} class object.
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
##' @export
bootstrap <- function(ylist, res){

  ## Setup
  numclust = res$numclust
  TT = length(ylist)
  dimdat = ncol(ylist %>% .[[1]])

  ## Conduct the E-step once to calculate responsibilities
  resp <- Estep(res$mn,
                res$sigma,
                res$prob,
                ylist = ylist,
                numclust = res$numclust,
                first_iter = TRUE)

  ## Do a coin-flip draw of membership
  drawslist = draw_membership(resp)

  ## Draw residuals
  residuals = get_residuals(ylist, res, drawslist)

  ## Pool all the residuals
  residuals_pooled = lapply(1:res$numclust, function(iclust){
    do.call(rbind, lapply(1:TT, function(tt) residuals[[tt]][[iclust]]))
  })

  ## Get the number of coin-flipped particles for each cluster
  ntklist = sapply(residuals, function(resids) resids %>% sapply(., nrow)) %>% t()


  ## Draw the new particles.
  new_ylist <- lapply(1:TT, function(tt){
    new_y_by_clust <- lapply(1:numclust, function(iclust){
      ntk = ntklist[tt, iclust]
      if(ntk == 0) return(rep(NA,dimdat) %>% rbind() %>% .[-1,])
      bootrows = sample(1:sum(ntklist[,iclust]), ntk, replace = TRUE)
      resampled_resids <- residuals_pooled[[iclust]][bootrows,,drop=FALSE]
      new_y <- resampled_resids %>% sweep(., 2, res$mn[tt,,iclust], "+")
      return(new_y)
    })
    return(new_y_by_clust)
  })
  new_ylist %<>% lapply(., function(a) do.call(rbind, a))

  return(list(ylist = new_ylist,
              X = res$X,
              ntklist = ntklist,
              drawslist = drawslist,
              residuals = residuals))
}



##' Drawing memberships by coinflip using the responsibilities.
##'
##' @importFrom magrittr %>%
##' @export
draw_membership <- function(resp){
  TT = length(resp)
  drawslist = lapply(1:TT, function(tt){
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
##' @inheritParams ylist
##' @inheritParams res
##'
##' @export
get_residuals <- function(ylist, res, drawslist){
  numclust = res$numclust
  mn = res$mn
  TT = length(ylist)
  residuals_by_cluster = lapply(1:TT, function(tt){
    y = ylist[[tt]]
    draws = drawslist[[tt]]
    resid_by_clust = lapply(1:numclust, function(iclust){
      this_clust_y = y[which(draws[,iclust]==1),,drop=FALSE]
      this_clust_residuals = sweep(this_clust_y, 2, mn[tt,,iclust])
    })
    names(resid_by_clust) = paste0("Clust", 1:res$numclust)
    return(resid_by_clust)
  })
  return(residuals_by_cluster)
}


##' Checks whether the dataset and \code{flowmix}-class object are
##' compatible. (Under construction! Can you think of anything else?)
##' @inheritParams ylist
##' @param res \code{flowmix} class object.
##' @importFrom assertthat assert_that
check_compatible <- function(ylist, res){
  assertthat::assert_that(res$TT == length(ylist))
  assertthat::assert_that(res$dimdat == ncol(ylist[[1]]))
}

##' Checks whether two datasets of \code{flowmix} class are exactly the same
##' size.
##' @param ylist1 One \code{flowmix} class object.
##' @param ylist2 Another \code{flowmix} class object.
##' @importFrom assertthat assert_that
check_if_same_size <- function(ylist1, ylist2){
  assertthat::assert_that(length(ylist1) == length(ylist2))
  assertthat::assert_that(all(ylist1$ntlist == ylist2$ntlist))
}
