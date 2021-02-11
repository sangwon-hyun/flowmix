##' Define the time folds cross-validation.
##'
##' @param blocksize Size (number of time indices) of one block.
##' @return List of fold indices.
##' @export
##'
make_cv_folds <- function(ylist, nfold, blocksize = 20){

  ## Make hour-long index list
  TT = length(ylist)
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

  ## plot(NA, xlim = c(0,TT), ylim=1:2)
  ## lapply(1:nfold, function(ifold){a = test.ii.list[[ifold]]; abline(v=a, col=ifold)})

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
                    ## The rest that is needed explicitly for flowmix()
                    ylist, countslist,
                    X, ...){

  ## Get the train/test data
  test.inds = unlist(folds[ifold])
  test.dat = ylist[test.inds]
  test.count = countslist[test.inds]
  test.X = X[test.inds,]
  train.dat = ylist[-test.inds]
  train.count = countslist[-test.inds]
  train.X = X[-test.inds,]

  ## Check whether this job has been done already.
  filename = make_cvscore_filename(ialpha, ibeta, ifold, irep, sim, isim)
  if(file.exists(file.path(destin, filename))){
    cat(filename, "already done", fill=TRUE)
    return(NULL)
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
                          ## The rest that is needed explicitly for flowmix_once()
                          ylist, countslist, X,
                          sim = FALSE, isim = 1,
                          ...){

  args = list(...)
  nrep = args$nrep
  for(irep in 1:nrep){

    ## Writing file
    filename = make_refit_filename(ialpha, ibeta, irep, sim, isim)
    if(file.exists(file.path(destin, filename))){
      cat(filename, "already done", fill=TRUE)
      next
    } else {

      ## Get the fitted results on the entire data
      args = list(...)
      args$ylist = ylist
      args$countslist = countslist
      args$X = X
      args$prob_lambda = prob_lambdas[ialpha]
      args$mean_lambda = mean_lambdas[ibeta]
      if("nrep" %in% names(args)) args = args[-which(names(args) %in% "nrep")] ## remove |nrep| prior to feeding

      ## Call the function.
      argn <- lapply(names(args), as.name)
      names(argn) <- names(args)
      call <- as.call(c(list(as.name("flowmix_once")), argn))
      res = eval(call, args)

      ## Save the results
      cat("Saving file here:", file.path(destin, filename), fill=TRUE)
      save(res, args, objective, file=file.path(destin, filename))
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
make_cvscore_filename <- function(ialpha, ibeta, ifold, irep,
                                  ## If simulations, then additional file names
                                  sim = FALSE, isim = 1){
  filename = paste0(ialpha, "-", ibeta, "-", ifold, "-", irep, "-cvscore.Rdata")
  if(sim){filename = paste0(isim, "-", filename)} ## Temporary
  return(filename)
}


##' Create file name (a string) for re-estimated models for the lambda values
##' indexed by \code{ialpha} \ibeta{ibeta}.
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
##' @param ylist data
##' @param X covariates
##' @param countslist multiplicity
##' @param destin Where to save the output.
##' @param max_mean_lambda Maximum value of regularization
##' @param nfold Number of cross-validation folds. Defaults to 5.
##' @param cv_gridsize Grid size for cross validation.
##' @param cv_fold_blocksize Number of time blocks to be used for cross-validation folds.
##' @param nrep Number of repetitions.


##' @param ... Additional arguments to flowmix().
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

  ## Define the CV folds (e.g. 5 big consecutive time blocks)
  ## folds = blockcv_make_folds(ylist = ylist, nfold = 5)
  folds = make_cv_folds(ylist = ylist, nfold = nfold, blocksize = 1)

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
           ## Save the file
           file = file.path(destin, 'meta.Rdata'))
      print(paste0("wrote meta data to ", file.path(destin, 'meta.Rdata')))
    }
  }

  ## Run the EM algorithm many times, for each value of (ialpha, ibeta, ifold, irep)
  start.time = Sys.time()
  parallel::mclapply(1:nrow(iimat), function(ii){
    printprogress(ii, nrow(iimat), "jobs (EM replicates), assigned on this computer", start.time = start.time)

    if(!refit){
      ialpha = iimat[ii,"ialpha"]
      ibeta = iimat[ii,"ibeta"]
      ifold = iimat[ii,"ifold"]
      irep = iimat[ii,"irep"]
      ## if(verbose) cat('(ialpha, ibeta, ifold, irep)=', c(ialpha, ibeta, ifold, irep), fill=TRUE)
    } else {
      ialpha = iimat[ii, "ialpha"]
      ibeta = iimat[ii, "ibeta"]
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
              verbose = FALSE)
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
                    verbose = FALSE)
    }
    return(NULL)
  }, mc.cores = mc.cores)
}
