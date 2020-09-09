##' Synopsis: BLOCKED cross-validation function goes here.


##' CV wrapper for flowmix().
##' @param nsplit Number of CV splits. Defaults to 5.
##' @param folds CV folds.
##' @param cv_gridsize Grid size of CV.
##' @param destin Destination directory.
##' @param mean_lambdas List of regularization parameters for mean model.
##' @param pie_lambdas List of regularization parameters for pie model.
##' @param iirange Which of the \code{cv_gridsize^2*nfold*nrep} jobs (e.g. of
##'   the form 5-3-2-5) to run. Defaults to \code{NULL}, in which case all jobs
##'   are run in order of these indices i.e. 1-1-1-1, then 1-1-1-2, and so on.
##' @param ... default arguments to flowmix().
##' @return List containing (1) the set of coefficients
##' @export
blockcv <- function(cl, folds, cv_gridsize, iirange=NULL,
                    ylist, countslist, X,
                    mean_lambdas,
                    pie_lambdas,
                    destin,
                    parallel=TRUE,
                    sim=FALSE,
                    isim=1,
                    ...){

  ## First, save the meta information.
  nfold = length(folds)
  args = list(...)

  ## |iimat| looks like this:
  ## (#, ialpha, ibeta, ifold, irep)
  ## 1, 1, 1, 1, 1
  ## 2, 1, 1, 1, 2
  ## 3, 1, 1, 1, 3
  ## ...
  ## 5625, 15, 15, 5, 5
  nfold = length(folds)
  nrep = args$nrep
  iimat = make_iimat(cv_gridsize, nfold, nrep)
  iimax = nrow(iimat)
  if(is.null(iirange)){
    iirange = 1:iimax
  }
  if(parallel){
    parallel::parLapplyLB(cl, iirange, function(ii){
      ialpha = iimat[ii,"ialpha"]
      ibeta = iimat[ii,"ibeta"]
      ifold = iimat[ii,"ifold"]
      irep = iimat[ii,"irep"]

      ## Temporary: If simulation
      if(sim){
      cat('(isim, ialpha, ibeta, ifold, irep)=', c(sim, ialpha, ibeta, ifold, irep), fill=TRUE)
      } else {
      cat('(ialpha, ibeta, ifold, irep)=', c(ialpha, ibeta, ifold, irep), fill=TRUE)
      }
      one_job(ialpha, ibeta, ifold, irep,
              folds, destin, mean_lambdas, pie_lambdas,
              ylist, countslist, X,
              sim, isim,## Temporary
              ...)
      cat(fill=TRUE)
    })
  } else {
    ## Mirror copy of whatever goes in the parallel=TRUE block directly above.
    lapply(iirange, function(ii){
      ialpha = iimat[ii,"ialpha"]
      ibeta = iimat[ii,"ibeta"]
      ifold = iimat[ii,"ifold"]
      irep = iimat[ii,"irep"]
      cat('iii=', c(ialpha, ibeta, ifold, irep), fill=TRUE)
      one_job(ialpha, ibeta, ifold, irep,
              folds, destin, mean_lambdas, pie_lambdas,
              ylist, countslist,
              sim, isim,## Temporary
              X, ...)
      cat(fill=TRUE)
    })
  }
}


##' CV wrapper for flowmix().
##'
##' @param nsplit Number of CV splits. Defaults to 5.
##' @param ... Other arguments to flowmix().
##'
##' @return List containing (1) the set of coefficients
##'
##' @export
blockcv_fitmodel <- function(cl, destin,
                             cv_gridsize,
                             mean_lambdas,
                             pie_lambdas,
                             ## Arguments to flowmix()
                             ylist,
                             countslist,
                             X,
                             ...){

  iimat = make_iimat_small(cv_gridsize)
  iimax = nrow(iimat)
  parallel::parLapplyLB(cl, 1:iimax, function(ii){
    ialpha = iimat[ii, "ialpha"]
    ibeta = iimat[ii, "ibeta"]
    print('iii='); print(c(ialpha, ibeta))
    one_job_refit(ialpha, ibeta, destin,
                  mean_lambdas, pie_lambdas,
                  ## Arguments to flowmix()
                  ylist, countslist, X,
                  ...)
  })
}

## ##' A wrapper for \code{blockcv_make_folds()} and
## ##' \code{blockcv_hourlong_make_folds()} that saves to and loads from an output
## ##' file. (Incomplete)
## ##'
## ##' @param filenaem
## ##'
## ##' @return
## ##'
## blockcv_make_folds_saveload <- function(ylist, nfold){
## }



##' Define the folds for nfold-CV.
blockcv_make_folds <- function(ylist, nfold, verbose=FALSE){

  if(verbose) print("Large consecutive time blocks used for CV (block type 1)")
  TT = length(ylist)
  endpoints = round(seq(from = 0, to = TT, length = nfold+1))
  inds = Map(function(a,b){(a+1):b},
             endpoints[-length(endpoints)],
             endpoints[-1])
  return(inds)
}

##' Define the folds for 1-hr-split CVs. (roughly, 1 hour = block size of 20).
blockcv_hourlong_make_folds <- function(ylist, nfold, verbose=FALSE, blocksize=20){

  if(verbose) print("Hour-long time blocks used for CV (block type 2)")

  ## Make hour-long index list
  TT = length(ylist)
  endpoints = round(seq(from = 0, to = TT + blocksize,
                        by = blocksize))
  inds = Map(function(a,b){
    if(a>=TT) return(NULL)
    (a+1):pmin(b,TT)
  },
             endpoints[-length(endpoints)],
             endpoints[-1])
  null.elt = sapply(inds, is.null)
  if(any(null.elt)){
    inds = inds[-which(null.elt)]
  }

  ## Further make these into five blocks of test indices.
  test.ii.list = lapply(1:nfold, function(ifold){
    which.test.inds = seq(from=ifold, to=length(inds), by=5)
    test.ii = unlist(inds[which.test.inds])
    return(test.ii)
  })

  ## plot(NA, xlim = c(0,TT), ylim=1:2)
  ## lapply(1:nfold, function(ifold){a = test.ii.list[[ifold]]; abline(v=a, col=ifold)})

  return(test.ii.list)
}


##' Helper function to run ONE job for blockCV, in ialpha, ibeta, ifold, irep.
##'
##' @param ialpha Index for alpha.
##' @param ibeta Index for beta.
##' @param ifold Index for CV folds.
##' @param irep Index for 1 through nrep.
##' @param folds CV folds.
##' @param destin Destination directory.
##' @param mean_lambdas List of regularization parameters for mean model.
##' @param pie_lambdas List of regularization parameters for pie model.
##' @param ylist Data.
##' @param countslist Counts or biomass.
##' @param X Covariates.
##' @param ... Rest of arguments for \code{flowmix_once()}.
##'
##' @return Nothing is returned. Instead, a file named "1-1-1-1-cvscore.Rdata"
##'   is saved in \code{destin}. (The indices here are ialpha-ibeta-ifold-irep).
##'
one_job <- function(ialpha, ibeta, ifold, irep, folds, destin,
                    mean_lambdas, pie_lambdas,
                    ## The rest that is needed explicitly for flowmix()
                    ylist, countslist,
                    sim, isim,## Temporary
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

  pie_lambda = pie_lambdas[ialpha]
  mean_lambda = mean_lambdas[ibeta]

  ## Run the algorithm (all this trouble because of |nrep|)
  args = list(...)
  args$ylist = train.dat
  args$countslist = train.count
  args$X = train.X
  args$mean_lambda = mean_lambda
  args$pie_lambda = pie_lambda
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

    ## Assign mn and pie
    pred = predict.flowmix(res.train, newx = test.X)
    stopifnot(all(pred$newpie >= 0))

    ## Evaluate on test data, by calculating objective (penalized likelihood)
    cvscore = objective(mu = pred$newmn,
                        pie = pred$newpie,
                        sigma = pred$sigma,
                        ylist = test.dat,
                        countslist = test.count,
                        pie_lambda = 0,
                        mean_lambda = 0,
                        alpha = res.train$alpha,
                        beta = res.train$beta)

    time_per_iter = res.train$time_per_iter
    final_iter = res.train$final.iter
    total_time = res.train$total_time
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
         pie_lambda,
         mean_lambdas,
         pie_lambdas,
         beta,
         alpha,
         objectives,
         ## Save the file
         file = file.path(destin, filename))
    return(NULL)

  }, error = function(err) {
    err$message = paste(err$message,
                        "\n(No file will be saved for lambdas (",
                        signif(pie_lambdas[ialpha],3), ", ", signif(mean_lambdas[ibeta],3),
                        ") whose indices are: ",
                        ialpha, "-", ibeta, "-", ifold, "-", irep,
                        " .)",sep="")
    cat(err$message, fill=TRUE)
    warning(err)})
}


##' Refit one job
one_job_refit <- function(ialpha, ibeta, destin,
                          mean_lambdas, pie_lambdas,
                          ## The rest that is needed explicitly for flowmix_once()
                          ylist, countslist, X,
                          sim = FALSE, isim=1,
                          ...){

  args = list(...)
  nrep = args$nrep
  for(irep in 1:nrep){

    ## Writing file
    filename = make_refit_filename(ialpha, ibeta, irep, sim, isim)
    if(file.exists(file.path(destin, filename))){
      cat(filename, "already done", fill=TRUE)
      return(NULL)
    } else {

      ## Get the fitted results on the entire data
      args = list(...)
      args$ylist = ylist
      args$countslist = countslist
      args$X = X
      args$pie_lambda = pie_lambdas[ialpha]
      args$mean_lambda = mean_lambdas[ibeta]
      if("nrep" %in% names(args)) args = args[-which(names(args) %in% "nrep")] ## remove |nrep| prior to feeding

      ## Call the function.
      argn <- lapply(names(args), as.name)
      names(argn) <- names(args)
      call <- as.call(c(list(as.name("flowmix_once")), argn))
      res = eval(call, args)

      ## Save the results
      print(res)
      cat("Saving file here:", file.path(destin, filename), fill=TRUE)
      save(res, args, objective, file=file.path(destin, filename))
    }
  }
}

##' Indices for the cross validation jobs.
##' The resulting iimat looks like this:
##' (#, ialpha, ibeta, ifold)
##' 1, 1, 1, 1
##' 2, 1, 1, 2
##' 3, 1, 1, 3
##' ...
##' 1125, 15, 15, 5
##'
##' @param cv_gridsize CV grid size.
##' @param nfold Number of of CV folds.
##'
##' @return Numeric matrix.
##'
make_iimat <- function(cv_gridsize, nfold, nrep){
  iimat = expand.grid(ialpha = 1:cv_gridsize,
                      ibeta = 1:cv_gridsize,
                      ifold = 1:nfold,
                      irep = 1:nrep)
  iimat = cbind(as.numeric(rownames(iimat)), iimat)
  return(iimat)
}

##' 2d indices for the cross validation jobs.
##' The resulting iimat looks like this:
##' (#, ialpha, ibeta)
##' 1, 1, 1
##' 2, 1, 2
##' 3, 1, 3
##' ...
##' 225, 15, 15
##'
##' @param cv_gridsize CV grid size.
make_iimat_small <- function(cv_gridsize){
  iimat = expand.grid(ialpha = 1:cv_gridsize,
                      ibeta = 1:cv_gridsize)
  iimat = cbind(as.numeric(rownames(iimat)), iimat)
  return(iimat)
}


make_cvscore_filename <- function(ialpha, ibeta, ifold, irep, sim, isim){
  filename = paste0(ialpha, "-", ibeta, "-", ifold, "-", irep, "-cvscore.Rdata")
  if(sim){filename = paste0(isim, "-", filename)} ## Temporary
  return(filename)
}

make_refit_filename <- function(ialpha, ibeta, irep, sim, isim){
  filename = paste0(ialpha, "-", ibeta, "-", irep, "-fit.Rdata")
  if(sim){filename = paste0(isim, "-", filename)} ## Temporary
  return(filename)
}
