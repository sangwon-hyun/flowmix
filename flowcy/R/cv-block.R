##' Synopsis: BLOCKED cross-validation function goes here.


##' CV wrapper for covarem().
##' @param nsplit Number of CV splits. Defaults to 5.
##' @param ... default arguments to covarem().
##' @return List containing (1) the set of coefficients
##' @export
blockcv <- function(cl, folds, cv_gridsize,
                    ylist, countslist, X,
                    mean_lambdas,
                    pie_lambdas,
                    destin,
                    parallel=TRUE,
                    ...){

  ## iimat looks like this:
  ## (#, ialpha, ibeta, ifold)
  ## 1, 1, 1, 1
  ## 2, 1, 1, 2
  ## 3, 1, 1, 3
  ## ...
  ## 1125, 15, 15, 5
  nfold = length(folds)
  iimat = make_iimat(cv_gridsize, nfold)
  iimax = nrow(iimat)
  if(parallel){
    parallel::parLapplyLB(cl,
    ## lapply(
        1:iimax, function(ii){
      ialpha = iimat[ii,"ialpha"]
      ibeta = iimat[ii,"ibeta"]
      ifold = iimat[ii,"ifold"]
      cat('(ialpha, ibeta, ifold)=', c(ialpha, ibeta, ifold), fill=TRUE)
      one_job(ialpha, ibeta, ifold, folds, destin,
              mean_lambdas, pie_lambdas,
              ylist, countslist, X, ...)
      cat(fill=TRUE)
    })
  } else {
    ## Mirror copy of whatever goes in the parallel=TRUE block directly above.
    lapply(1:iimax, function(ii){
      ialpha = iimat[ii,"ialpha"]
      ibeta = iimat[ii,"ibeta"]
      ifold = iimat[ii,"ifold"]
      cat('iii=', c(ialpha, ibeta, ifold), fill=TRUE)
      one_job(ialpha, ibeta, ifold, folds, destin,
              mean_lambdas, pie_lambdas,
              ylist, countslist, X, ...)
      cat(fill=TRUE)
    })
  }
}


##' CV wrapper for covarem().
##' @param nsplit Number of CV splits. Defaults to 5.
##' @param ... default arguments to covarem().
##' @return List containing (1) the set of coefficients
##' @export
blockcv_fitmodel <- function(cl, destin,
                             ylist, countslist, X,
                             mean_lambdas,
                             pie_lambdas,
                             ...){

  args = list(...)
  iimat = make_iimat_small(args$cv_gridsize)
  iimax = nrow(iimat)
  parallel::parLapplyLB(cl,
  ## 1:iimax, function(ii){
  ## lapply(
      1:iimax, function(ii){
    ialpha = iimat[ii, "ialpha"]
    ibeta = iimat[ii, "ibeta"]
    print('iii=')
    print(c(ialpha, ibeta))
    one_job_refit(ialpha, ibeta, destin,
                  mean_lambdas, pie_lambdas,
                  ylist, countslist, X,
                  ...)
  })

}




##' Define the folds for nfold-CV.
blockcv_make_folds <- function(ylist, nfold, verbose=FALSE){

  if(verbose) print("Large consecutive time blocks used for CV (block type 1)")
  TT = length(ylist)
  endpoints = round(seq(from = 1, to = TT, length = nfold+1))
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
  endpoints = round(seq(from = 0, to = TT + blocksize, by = blocksize))
  inds = Map(function(a,b){(a+1):b},
             endpoints[-length(endpoints)],
             endpoints[-1])

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


##' Helper function to run ONE job for blockCV
one_job <- function(ialpha, ibeta, ifold, folds, destin,
                    mean_lambdas, pie_lambdas,
                    ## The rest that is needed explicitly for covarem()
                    ylist, countslist, X, ...){

  ## Get the train/test data
  test.inds = unlist(folds[ifold])
  test.dat = ylist[test.inds]
  test.count = countslist[test.inds]
  test.X = X[test.inds,]
  train.dat = ylist[-test.inds]
  train.count = countslist[-test.inds]
  train.X = X[-test.inds,]


  ## Check whether this version has been done already.
  filename = paste0(ialpha, "-", ibeta, "-", ifold, "-cvscore.Rdata")
  if(file.exists(file.path(destin, filename))){
    cat("ialpha, ibeta ", ialpha, ibeta, "are already done.", fill=TRUE)
    return(NULL)
  }
  mean_lambda = mean_lambdas[ibeta]
  pie_lambda = pie_lambdas[ialpha]

  ## tryCatch({

    ## Run algorithm on training data,
    res.train = covarem(ylist = train.dat,
                        countslist = train.count,
                        X = train.X,
                        ## refit = FALSE, ## Is this necessary? Think about it for a sec.
                        mean_lambda = mean_lambda,
                        pie_lambda = pie_lambda,
                        ## verbose=TRUE, ## TEMPORARY
                        ## plot = TRUE,
                        ## plotdir = paste0("~/Desktop/blockcv-test-figures/1-6"),
                        ## filepath = file.path("~/Desktop/blockcv-test-figures/1-6", Sys.time()),
                        ...)
    ## assert_that(!refit)

    ## Assign mn and pie
    pred = predict.covarem(res.train, newx = test.X)
    stopifnot(all(pred$newpie >= 0))

    ## Evaluate on test data, by calculating objective (penalized likelihood)
    cvscore = objective(mu = pred$newmn,
                        pie = pred$newpie,
                        sigma = pred$sigma,
                        ylist = test.dat,
                        countslist = NULL,
                        pie_lambda = 0,
                        mean_lambda = 0,
                        alpha = res.train$alpha,
                        beta = res.train$beta)

    time_per_iter = res.train$time_per_iter
    final_iter = res.train$final.iter
    total_time = res.train$total_time
    beta = res.train$beta
    alpha = res.train$alpha

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
         ## Save the file
         file = file.path(destin, filename))

    return(NULL)

  ## }, error = function(err) {
  ##   err$message = paste(err$message,
  ##                       "\n(No file will be saved for the lambdas ",
  ##                       pie_lambdas[ialpha], ",", mean_lambdas[ibeta],
  ##                       "whose indices are", ialpha, "-", ibeta, "-", ifold,
  ##                       " .)",sep="")
  ##   cat(err$message, fill=TRUE)
  ##   warning(err)})
}



#############################
### Helper to run one job ###
#############################

##' Refit one job
one_job_refit <- function(ialpha, ibeta, destin,
                          mean_lambdas, pie_lambdas,
                          ## The rest that is needed explicitly for covarem()
                          ylist, countslist, X,
                          ...){

  ## Get the fitted results on the entire data
  res = covarem(ylist = ylist, countslist = countslist, X = X,
                mean_lambda = mean_lambdas[ibeta],
                pie_lambda = pie_lambdas[ialpha],
                ...)

  ## Save the results
  filename = paste0(ialpha, "-", ibeta, "-fit.Rdata")
  save(res, file=file.path(destin, filename))
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
make_iimat <- function(cv_gridsize, nfold){
  iimat = expand.grid(ialpha = 1:cv_gridsize,
                      ibeta = 1:cv_gridsize,
                      ifold = 1:nfold)
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



##' Aggregate results from blocked CV.
##'
##' @param destin destination folder for output.
##' @param cv_gridsize Grid size for cross validation.
##' @param nfold Number of folds.
##'
##' @export
blockcv_aggregate <- function(destin, cv_gridsize, nfold){
  cvscore.array = array(0, dim=c(cv_gridsize, cv_gridsize, nfold))
  cvscore.mat = matrix(0, nrow = cv_gridsize, ncol = cv_gridsize)
  for(ialpha in 1:cv_gridsize){
    for(ibeta in 1:cv_gridsize){
      cvscores <- sapply(1:nfold, function(igrid){
        filename = paste0(ialpha, "-", ibeta, "-", igrid, "-cvscore.Rdata")
        print(filename)
        tryCatch({
          load(file.path(destin, filename))
          return(cvscore)
        }, error = function(e){
          return(NA)
        })
      })
      cvscore.mat[ialpha, ibeta] = mean(cvscores)
      cvscore.array[ialpha, ibeta, ] = cvscores
    }
  }
  return(list(cvscore.array = cvscore.array,
              cvscore.mat = cvscore.mat))
}

