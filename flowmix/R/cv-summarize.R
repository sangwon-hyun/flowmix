##' Main function for summarizing the cross-validation results.
##'
##' @inheritParams cv.flowmix
##' @param filename File name to save to.
##' @param save If TRUE, save to \code{file.path(destin, filename)}.
##'
##' @return List containing various outcomes from the cross-validation, such as
##'   \code{bestres} which is the \code{flowmix} class object of the overall
##'   best model chosen from the cross-validation; \code{cvscoremat} containing
##'   a 2d matrix of CV scores from all pairs of lambdas; \code{bestreslist}
##'   contains all the best models (out of \code{nrep} EM replications)q from the
##'   each pair of lambda values. If \code{isTRUE(save)}, nothing is returned.
##'
##' @export
cv_summary <- function(destin = ".",
                       save = FALSE,
                       filename = "summary.RDS"
                       ){

  ####################
  ## Load data #######
  ####################
  load(file = file.path(destin, 'meta.Rdata'), verbose = FALSE)

  ## This loads all the necessary things: nrep, nfold, cv_gridsize
  stopifnot(exists("nrep"))
  stopifnot(exists("nfold"))
  stopifnot(exists("cv_gridsize"))

  ##########################
  ## Get the CV results. ###
  ##########################
  a = cv_aggregate(destin)
  cvscore.mat = a$cvscore.mat
  min.inds = a$min.inds

  ## Get the refit flowmix results
  bestreslist = cv_aggregate_res(destin = destin)
  bestres = bestreslist[[paste0(min.inds[1] , "-", min.inds[2])]]
  if(is.null(bestres)){
    stop(paste0("The model with lambda indices (",
                min.inds[1], ",", min.inds[2], ") is not available."))
  }

  if(is.null(colnames(bestres$X))){
    colnames(bestres$X) = 1:ncol(bestres$X)
  }

  ########################
  ## Get coefficients ####
  ########################
  betalist =  lapply(1:bestres$numclust, function(iclust){
    ## Get all betas
    rownames(bestres$beta[[iclust]])[-1] = colnames(bestres$X)
    cf = bestres$beta[[iclust]][-1,, drop=FALSE]
    ## TODO: TRY dplyr here:

    ## Remove the rows that are all zero
    all.zero.rows = which(apply(cf, 1, function(myrow)all(myrow == 0)))
    if(length(all.zero.rows) > 0){
      cf = cf[-all.zero.rows,, drop=FALSE]
    }
    round(Matrix::Matrix(cf, sparse=TRUE),3)
  })
  names(betalist) = paste0("Beta matrix, cluster ", 1:bestres$numclust)
  pretty.betas = betalist
  colnames(bestres$alpha)[-1 ] = colnames(bestres$X)
  alpha = t(bestres$alpha)
  alpha[which(abs(alpha) < 1E-5)] = 0
  pretty.alphas = round(Matrix::Matrix(alpha, sparse=TRUE),3)

  ######################
  ## Get the sigmas ####
  ######################
  if(bestres$dimdat == 1){
    pretty.sigmas = sqrt(bestres$sigma[,1,])
    names(pretty.sigmas) = paste0("Cluster ", 1:bestres$numclust)
  } else {
    sigmas = lapply(1:bestres$numclust, function(iclust){
      diag(bestres$sigma[iclust,,])
    })
    names(sigmas) = paste0("Cluster ", 1:bestres$numclust)
    pretty.sigmas = lapply(sigmas, sqrt)
  }

  out = list(bestres = bestres,
             cvscore.mat = cvscore.mat,
             min.inds = min.inds,
             ## Pretty formatted data
             pretty.alphas = pretty.alphas,
             pretty.betas = pretty.betas,
             pretty.sigmas = pretty.sigmas,
             ## List of all best models for all lambda pairs.
             bestreslist = bestreslist,
             destin = destin)

  if(save){ saveRDS(out, file=file.path(destin, filename)); invisible(NULL) }
  return(out)
}

##' Aggregate CV scores from the results, saved in \code{destin}.
##'
##' @param destin Directory with cross-validation output.
##' @param sim Simulation or not?
##' @param isim Simulation number.
##'
##' @export
cv_aggregate <- function(destin,
                         sim = FALSE,
                         isim = 1){

  ## ## Read the meta data (for |nfold|, |cv_gridsize|, |nrep|, |prob_lambdas|,
  ## ## |mean_lambdas|)
  load(file = file.path(destin, 'meta.Rdata'), verbose = FALSE)

  ## This loads all the necessary things; just double-checking.
  stopifnot(exists("nrep"))
  stopifnot(exists("nfold"))
  stopifnot(exists("cv_gridsize"))
  stopifnot(exists(c("prob_lambdas", "mean_lambdas")))

  ## Purely for back-compatability (retire soon)
  ## if(exists("pie_lambdas")) prob_lambdas = pie_lambdas

  ## Aggregate the results
  cvscore.array = array(NA, dim = c(cv_gridsize, cv_gridsize, nfold, nrep))
  cvscore.mat = matrix(NA, nrow = cv_gridsize, ncol = cv_gridsize)
  for(ialpha in 1:cv_gridsize){
    for(ibeta in 1:cv_gridsize){
      obj = matrix(NA, nrow=nfold, ncol=nrep)
      for(ifold in 1:nfold){
        for(irep in 1:nrep){
          filename = make_cvscore_filename(ialpha, ibeta, ifold, irep, sim, isim)
          tryCatch({
            load(file.path(destin, filename), verbose = FALSE)


            cvscore.array[ialpha, ibeta, ifold, irep] = cvscore
            obj[ifold, irep] = objectives[length(objectives)]
          }, error = function(e){})
        }
      }

      ## Pick out the CV scores with the *best* (lowest) objective value
      cvscores = cvscore.array[ialpha, ibeta,,]
      best.models = apply(obj, 1, function(myrow){
        ind = which(myrow == min(myrow, na.rm=TRUE))
        if(length(ind)>1) ind = ind[1]  ## Just choose one, if there is a tie.
        return(ind)
      })
      final.cvscores = sapply(1:nfold, function(ifold){
        cvscores[ifold, best.models[ifold]]
      })
      cvscore.mat[ialpha, ibeta] = mean(final.cvscores)
    }
  }

  ## Clean a bit
  cvscore.mat[which(is.nan(cvscore.mat), arr.ind=TRUE)] = NA

  ## ## Read the meta data (for |nfold|, |cv_gridsize|, |nrep|)
  rownames(cvscore.mat) = signif(prob_lambdas,3)
  colnames(cvscore.mat) = signif(mean_lambdas,3)


  ## Find the minimum
  mat = cvscore.mat
  min.inds = which(mat == min(mat, na.rm = TRUE), arr.ind = TRUE)

  ## Return the results
  out = list(cvscore.array = cvscore.array,
              cvscore.mat = cvscore.mat,
              mean_lambdas = mean_lambdas,
              prob_lambdas = prob_lambdas,
              min.inds = min.inds)
  return(out)
}




##' Helper to aggregate parallelized CV results and obtain the |res| object, all
##' saved in |destin|.
##'
##' @inheritParams cv_aggregate
##'
##' @return List containing, for every (ialpha, ibeta), the "best" estimated
##'   model out of the |nrep| replicates (best in the sense that it had the best
##'   likelihood value out of the |nrep| replicates.)
cv_aggregate_res <- function(destin,
                             ## Is this a simulation or not? (soon to be outdated)
                             sim = FALSE,
                             isim = NULL
                             ){

  load(file.path(destin, "meta.Rdata"))

  ## df.mat = matrix(NA, ncol=cv_gridsize, nrow=cv_gridsize)
  res.list = list()
  for(ialpha in 1:cv_gridsize){
    for(ibeta in 1:cv_gridsize){

      ## Objective values, over nrep
      obj = rep(NA, nrep)
      ## df = rep(NA, nrep) #
      res.list.inner = list()
      for(irep in 1:nrep){
        ## filename = paste0(ialpha, "-", ibeta, "-", irep, "-fit.Rdata")
        filename = make_refit_filename(ialpha, ibeta, irep, sim, isim)
        tryCatch({
          ## Load fitted result
          load(file.path(destin, filename))

          ## Calculate DF
          res.list.inner[[irep]] = res

          ## Also calculate objective function
          obj[irep] = res$objectives[length(res$objectives)]

        }, error = function(e){})
      }

      ## Calculate the df of the best model
      if(!all(is.na(obj))){
        res.list[[paste0(ialpha, "-", ibeta)]] = res.list.inner[[which.min(obj)]] ## which.min?
      }
    }
  }
  return(res.list)
}




