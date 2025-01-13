##' Define the time folds for cross-validation.
##'
##' @param blocksize Size (number of time indices) of one block.
##' @param ylist Data.
##' @param nfold Number of CV folds.
##' @param blocksize Size of time blocks for forming groups of time points
##'   (cross-validation folds).
##' @param TT Number of time points in total.
##'
##' @return List of fold indices.
##' @export
##'
make_cv_folds <- function(ylist=NULL, nfold, blocksize = 20, TT=NULL){

  ## Make hour-long index list
  if(is.null(TT)) TT = length(ylist)

  endpoints = round(seq(from = 0, to = TT + blocksize,
                        by = blocksize))
  inds = Map(function(a,b){
    if(a >= TT) return(NULL)
    ## (a+1):pmin(b,TT))
    return(seq(a+1, pmin(b,TT)))
  }, endpoints[-length(endpoints)],  endpoints[-1])
  null.elt = sapply(inds, is.null)
  if(any(null.elt)){
    inds = inds[-which(null.elt)]
  }

  ## Further make these into five blocks of test indices.
  test.ii.list = lapply(1:nfold, function(ifold){
    which.test.inds = seq(from = ifold, to = length(inds), by = nfold)
    test.ii = unlist(inds[which.test.inds])
    return(test.ii)
  })

  ## Useful plotting code showing the plots.
  if(FALSE){
    plot(NA, xlim = c(0,TT), ylim=1:2)
    lapply(1:nfold, function(ifold){a = test.ii.list[[ifold]]; abline(v=a, col=ifold)})
  }

  return(test.ii.list)
}

##' Using subset of indices \code{subsample_inds}, do what \code{make_cv_folds}
##' does, but using original memberships from data of size \code{orig_TT}. The
##' new memberships are a subset of those original memberships.
##'
##' @param nfold Number of CV folds.
##' @param blocksize Size (number of time indices) of one block.
##' @param orig_TT Original data size (for forming original memberships).
##' @param subsample_inds Indices, out of \code{1:orig_TT}.
##'
##' @export
make_cv_folds_subsample_with_original_membership <- function(nfold, blocksize, orig_TT, subsample_inds){

  ## ## Sample settings
  ## la('flowmix')
  ## orig_TT = 298
  ## blocksize = 20
  ## nfold = 5
  ## subsample_size = 100
  ## ## End of sample settings

  ## Basic checks
  stopifnot(all(subsample_inds %in% 1:orig_TT))

  ## Make original folds
  orig_folds = make_cv_folds(TT = orig_TT,##dat$ylist,
                             nfold = nfold,
                             blocksize = blocksize)
  membership = rep(NA, orig_TT)
  for(ii in 1:nfold){   membership[orig_folds[[ii]]] = ii }

  ## Test indices in the subsamples
  ## subsample_inds = sample(1:orig_TT, subsample_size, replace = FALSE) %>% sort()
  subsample_membership = membership[subsample_inds]
  test.ii.list <- lapply(1:nfold, function(ifold){
    which(subsample_membership == ifold)
  })

  ## ## Test indices in the full data
  ## test.ii.list <- lapply(1:nfold, function(ifold){
  ##   subsample_inds[test.ii.list[[ifold]]]
  ## })

  ## Basic checks
  stopifnot(length(test.ii.list) == nfold)
  stopifnot(all((test.ii.list %>% unlist() %>% sort()) == (1:length(subsample_inds))))

  return(test.ii.list)
}



##' Helper function to run ONE job for CV, in ialpha, ibeta, ifold, irep.
##'
##' @param ialpha Index for alpha.
##' @param ibeta Index for beta.
##' @param ifold Index for CV folds.
##' @param irep Index for 1 through nrep.
##' @param folds CV folds.
##' @param destin Destination directory.
##' @param mean_lambdas List of regularization parameters for mean model.
##' @param prob_lambdas List of regularization parameters for prob model.
##' @param ylist Data.
##' @param countslist Counts or biomass.
##' @param X Covariates.
##' @param sim Set to \code{TRUE} if this is a simulation; the file name
##'   containing results for, say `isim=3`, has "3-" appended to the beginning.
##' @param isim Simulation number.
##' @param seedtab A table containing seeds (7 columns)
##' @param ... Rest of arguments for \code{flowmix_once()}.
##'
##' @return Nothing is returned. Instead, a file named "1-1-1-1-cvscore.Rdata"
##'   is saved in \code{destin}. (The indices here are ialpha-ibeta-ifold-irep).
##'
##' If this is a simulation, the file name containing results for, say `isim=3`,
##'   has "3-" appended to the beginning.
##'
##' @export
one_job <- function(ialpha, ibeta, ifold, irep, folds, destin,
                    mean_lambdas, prob_lambdas,
                    ## For simulations
                    sim = FALSE, isim = 1,
                    seedtab = NULL,
                    ## The rest that is needed explicitly for flowmix()
                    ylist, countslist,
                    X,
                    ...){

  ## Solving warning from devtools::check()
  seed1 = seed2 = seed3 = seed4 = seed5 = seed6 = seed7 = NULL

  ## Get the train/test data
  test.inds = unlist(folds[ifold])
  # test.dat = ylist[test.inds]
  # test.count = countslist[test.inds]
  # test.X = X[test.inds,,drop=FALSE]
  # train.dat = ylist[-test.inds]
  # train.count = countslist[-test.inds]
  # train.X = X[-test.inds,, drop=FALSE]

  # Three cases: 
  # 1. names(ylist) & others are arbitrary
  # 2. names(ylist) & others are equal to as.character(1:length(ylist)) (this is to 
  #   allow cross-validation using a subset of the original data)
  # 3. at least one of the names() or rownames() is NULL

  # In the first and third case, we use test.inds directly to subset the data
  # In the second case, we assume that ylist, countslist, and X are subsets of the 
  #   original data but `test.inds` refers to indices from the original (entire) 
  #   dataset. names(ylist) and the like store 1:length(orig_ylist)
  #   so that subsetting does not change the indices of the data points.

  # Rough solution for now
  if(identical(names(ylist), 1:length(ylist)) & 
     identical(names(countslist), 1:length(countslist)) & 
     identical(rownames(X), 1:nrow(X))) # case 2
  {
    test.inds <- as.character(test.inds)
    train.inds <- setdiff(names(ylist), test.inds) 

    train.dat = ylist[train.inds]
    train.count = countslist[train.inds]
    train.X = X[train.inds, , drop = FALSE]
  } else { # cases 1 & 3
    train.dat = ylist[-test.inds]
    train.count = countslist[-test.inds]
    train.X = X[-test.inds,, drop=FALSE]
  }
  
  # cases 1, 2, and 3
  test.dat = ylist[test.inds]
  test.count = countslist[test.inds]
  test.X = X[test.inds,,drop=FALSE]

  ## Check whether this job has been done already.
  filename = make_cvscore_filename(ialpha, ibeta, ifold, irep, sim, isim)
  if(file.exists(file.path(destin, filename))){
    cat(filename, "already done", fill=TRUE)
    return(NULL)
  }


  ## Get the seed ready
  if(!is.null(seedtab)){
    seed = seedtab %>%
      dplyr::filter(ialpha == !!ialpha,
                    ibeta == !!ibeta,
                    ifold == !!ifold,
                    irep == !!irep) %>%
      dplyr::select(seed1, seed2, seed3, seed4, seed5, seed6, seed7) %>% unlist() %>% as.integer()
  } else {
    seed = NULL
  }


  prob_lambda = prob_lambdas[ialpha]
  mean_lambda = mean_lambdas[ibeta]

  ## Run the algorithm (all this trouble because of |nrep|)
  args = list(...)
  args$ylist = train.dat
  args$countslist = train.count
  args$X = train.X
  args$mean_lambda = mean_lambda
  args$prob_lambda = prob_lambda
  args$seed = seed
  if("nrep" %in% names(args)){
    args = args[-which(names(args) %in% "nrep")] ## remove |nrep| prior to feeding to flowmix_once().
  }

  tryCatch({

    ## New and better do.call() statement:
    ## res.train = do.call(flowmix_once, args) ## Old
    argn <- lapply(names(args), as.name)
    names(argn) <- names(args)
    call <- as.call(c(list(as.name("flowmix_once")), argn))
    res.train = eval(call, args)

    ## Assign mn and prob
    pred = predict.flowmix(res.train, newx = test.X)
    stopifnot(all(pred$prob >= 0))

    ## Evaluate on test data, by calculating objective (penalized likelihood)
    cvscore = objective(mu = pred$mn,
                        prob = pred$prob,
                        sigma = pred$sigma,
                        ylist = test.dat,
                        countslist = test.count,
                        prob_lambda = 0,
                        mean_lambda = 0,
                        alpha = res.train$alpha,
                        beta = res.train$beta)

    ## Store (temporarily) the run times
    time_per_iter = res.train$time_per_iter
    final_iter = res.train$final.iter
    total_time = res.train$total_time

    ## Store the results.
    beta = res.train$beta
    alpha = res.train$alpha
    objectives = res.train$objectives

    ## Save the CV results
    save(cvscore,
         ## Time
         time_per_iter,
         final_iter,
         total_time,
         ## Results
         mean_lambda,
         prob_lambda,
         mean_lambdas,
         prob_lambdas,
         beta,
         alpha,
         objectives,
         ## Save the file
         file = file.path(destin, filename))
    return(NULL)

  }, error = function(err) {
    err$message = paste(err$message,
                        "\n(No file will be saved for lambdas (",
                        signif(prob_lambdas[ialpha],3), ", ", signif(mean_lambdas[ibeta],3),
                        ") whose indices are: ",
                        ialpha, "-", ibeta, "-", ifold, "-", irep,
                        " .)",sep="")
    cat(err$message, fill=TRUE)
    warning(err)})
}


##' Refit model for one pair of regularization parameter values. Saves to
##' \code{nrep} files named like "1-4-3-fit.Rdata", for
##' "(ialpha)-(ibeta)-(irep)-fit.Rdata".
##'
##' (Note, \code{nrep} is not an input to this function.)
##'
##' @inheritParams one_job
##'
##' @export
one_job_refit <- function(ialpha, ibeta, destin,
                          mean_lambdas, prob_lambdas,
                          seedtab = NULL,
                          ## The rest that is needed explicitly for flowmix_once()
                          ylist, countslist, X,
                          sim = FALSE, isim = 1,
                          ...){

  ## Solving warning from devtools::check()
  seed1 = seed2 = seed3 = seed4 = seed5 = seed6 = seed7 = NULL

  args = list(...)
  nrep = args$nrep
  for(irep in 1:nrep){

    ## Writing file
    filename = make_refit_filename(ialpha, ibeta, irep, sim, isim)
    if(file.exists(file.path(destin, filename))){
      cat(filename, "already done", fill=TRUE)
      next
    } else {

     ## Get the seed ready
     if(!is.null(seedtab)){
       ifold = 0
       seed = seedtab %>%
         dplyr::filter(ialpha == !!ialpha,
                       ibeta == !!ibeta,
                       ifold == !!ifold,
                       irep == !!irep) %>%
         dplyr::select(seed1, seed2, seed3, seed4, seed5, seed6, seed7) %>% unlist() %>% as.integer()
     } else {
       seed = NULL
     }

      ## Get the fitted results on the entire data
      args = list(...)
      args$ylist = ylist
      args$countslist = countslist
      args$X = X
      args$prob_lambda = prob_lambdas[ialpha]
      args$mean_lambda = mean_lambdas[ibeta]
      args$seed = seed
      if("nrep" %in% names(args)) args = args[-which(names(args) %in% "nrep")] ## remove |nrep| prior to feeding

      ## Call the function.
      argn <- lapply(names(args), as.name)
      names(argn) <- names(args)
      call <- as.call(c(list(as.name("flowmix_once")), argn))
      res = eval(call, args)

      ## Save the results
      cat("Saving file here:", file.path(destin, filename), fill=TRUE)
      save(res, file=file.path(destin, filename))
    }
  }
}

##' Indices for the cross validation jobs.
##'
##' The resulting iimat looks like this:
##'
##' ind ialpha ibeta ifold irep
##'  55      6     1     2    1
##'  56      7     1     2    1
##'  57      1     2     2    1
##'  58      2     2     2    1
##'  59      3     2     2    1
##'  60      4     2     2    1
##' @param cv_gridsize CV grid size.
##' @param nfold Number of of CV folds.
##' @param nrep Number of repetitions.
##'
##' @return Integer matrix.
##'
##' @export
make_iimat <- function(cv_gridsize, nfold, nrep){
  iimat = expand.grid(ialpha = 1:cv_gridsize,
                      ibeta = 1:cv_gridsize,
                      ifold = 1:nfold,
                      irep = 1:nrep)
  iimat = cbind(ind = as.numeric(rownames(iimat)), iimat)
  return(iimat)
}

##' 2d indices for the cross validation jobs.
##'
##' The resulting iimat looks like this:
##' (#, ialpha, ibeta)
##' 1, 1, 1
##' 2, 1, 2
##' 3, 1, 3
##' ...
##' 225, 15, 15
##'
##' @param cv_gridsize CV grid size.
##'
##' @return Integer matrix.
##'
##' @export
make_iimat_small <- function(cv_gridsize){
  iimat = expand.grid(ialpha = 1:cv_gridsize,
                      ibeta = 1:cv_gridsize)
  iimat = cbind(as.numeric(rownames(iimat)), iimat)
  return(iimat)
}


##' Create file name (a string) for cross-validation results.
##' @noRd
make_cvscore_filename <- function(ialpha, ibeta, ifold, irep,
                                  ## If simulations, then additional file names
                                  sim = FALSE, isim = 1){
  filename = paste0(ialpha, "-", ibeta, "-", ifold, "-", irep, "-cvscore.Rdata")
  if(sim){filename = paste0(isim, "-", filename)} ## Temporary
  return(filename)
}


##' Create file name (a string) for re-estimated models for the lambda values
##' indexed by \code{ialpha} and \code{ibeta}.
##' @noRd
make_refit_filename <- function(ialpha, ibeta, irep,
                                  ## If simulations, then additional file names
                                sim = FALSE, isim = 1){
  filename = paste0(ialpha, "-", ibeta, "-", irep, "-fit.Rdata")
  if(sim){filename = paste0(isim, "-", filename)} ## Temporary
  return(filename)
}

##' Cross-validation wrapper for flowmix(). Saves results to separate files in
##' \code{destin}.
##'
##' @param destin Where to save the output.
##' @param nfold Number of cross-validation folds. Defaults to 5.
##' @param nrep Number of repetitions.
##' @param save_meta If TRUE, save meta data.
##' @param mean_lambdas Regularization parameters for betas.
##' @param prob_lambdas Regularization parameters for alphas.
##' @param folds Manually provide CV folds (list of time points of data to use
##'   as CV folds). Defaults to NULL.
##' @param mc.cores Use this many CPU cores.
##' @param blocksize Contiguous time blocks from which to form CV time folds.
##' @param refit If TRUE, estimate the model on the full data, for each pair of
##'   regularization parameters.
##' @param iimat Indices for the cross validation jobs -- made using
##'   \code{make_iimat()}.
##' @param maxdev Maximum deviation of cluster means across time.
##' @param seedtab A table containing seeds (7 columns)
##' @param verbose TRUE for loudness.

##' @param flatX_thresh Threshold for detecting if any covariates are flat (low
##'   variance). These flat coefficients will have be set to zero and excluded
##'   from estimation altogether.
##'
##' @param ... Additional arguments to flowmix().
##' @inheritParams flowmix_once
##'
##' @return No return.
##'
##' @export
cv.flowmix <- function(
                       ## Data
                       ylist,
                       countslist,
                       X,
                       ## Define the locations to save the CV.
                       destin = ".",
                       ## Regularization parameter values
                       mean_lambdas,
                       prob_lambdas,
                       iimat = NULL,
                       ## Other settings
                       maxdev,
                       numclust,
                       nfold,
                       nrep,
                       verbose = FALSE,
                       refit = FALSE,
                       save_meta = FALSE,
                       mc.cores = 1,
                       blocksize,
                       folds = NULL,
                       seedtab = NULL,
                       flatX_thresh = 1E-5,
                       ...){

  ## Basic checks
  stopifnot(length(prob_lambdas) == length(mean_lambdas))
  cv_gridsize = length(mean_lambdas)

  ## There's an option to input one's own iimat matrix.
  if(is.null(iimat)){
    ## Make an index of all jobs
    if(!refit) iimat = make_iimat(cv_gridsize, nfold, nrep)
    if(refit) iimat = make_iimat_small(cv_gridsize)
  }

  ## Define the CV folds
  ## folds = make_cv_folds(ylist = ylist, nfold = nfold, blocksize = 1)
  if(is.null(folds)){
    folds = make_cv_folds(ylist = ylist, nfold = nfold, blocksize = blocksize)
  } else {
    stopifnot(length(folds) == nfold)
  }

  ## Save meta information, once.
  if(save_meta){
    if(!refit){
      save(folds,
           nfold,
           nrep, ## Added recently
           cv_gridsize,
           mean_lambdas,
           prob_lambdas,
           ylist, countslist, X,
           ## bins, ## Not available anymore
           ## sparse_countslist, ## Not available anymore
           ## Save the file
           file = file.path(destin, 'meta.Rdata'))
      print(paste0("wrote meta data to ", file.path(destin, 'meta.Rdata')))
    }
  }

  ## Run the EM algorithm many times, for each value of (ialpha, ibeta, ifold, irep)
  start.time = Sys.time()
  parallel::mclapply(1:nrow(iimat), function(ii){
    print_progress(ii, nrow(iimat), "Jobs (EM replicates) assigned on this computer", start.time = start.time)

    if(!refit){
      ialpha = iimat[ii,"ialpha"]
      ibeta = iimat[ii,"ibeta"]
      ifold = iimat[ii,"ifold"]
      irep = iimat[ii,"irep"]
      ## if(verbose) cat('(ialpha, ibeta, ifold, irep)=', c(ialpha, ibeta, ifold, irep), fill=TRUE)
    } else {
      ialpha = iimat[ii, "ialpha"]
      ibeta = iimat[ii, "ibeta"]
      ifold = 0
    }

    if(!refit){
      ## Add noise to X, if applicable
      one_job(ialpha = ialpha,
              ibeta = ibeta,
              ifold = ifold,
              irep = irep,
              folds = folds,
              destin = destin,
              mean_lambdas = mean_lambdas,
              prob_lambdas = prob_lambdas,
              ## Arguments for flowmix()
              ylist = ylist, countslist = countslist, X = X,
              ## Additional arguments for covarem(), for ellipsis.
              numclust = numclust,
              maxdev = maxdev,
              verbose = FALSE,
              seedtab = seedtab,
              flatX_thresh = flatX_thresh)
    } else {
      one_job_refit(ialpha = ialpha,
                    ibeta = ibeta,
                    destin = destin,
                    mean_lambdas = mean_lambdas,
                    prob_lambdas = prob_lambdas,
                    ## Arguments to flowmix()
                    ylist = ylist, countslist = countslist, X = X,
                    ## Additional arguments for covarem(), for ellipsis.
                    numclust = numclust,
                    maxdev = maxdev,
                    nrep = nrep,
                    verbose = FALSE,
                    seedtab = seedtab,
                    flatX_thresh = flatX_thresh
                    )
    }
    return(NULL)
  }, mc.cores = mc.cores)
}

