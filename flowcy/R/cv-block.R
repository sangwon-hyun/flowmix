##' Synopsis: BLOCKED cross-validation function goes here.


##' CV wrapper for covarem().
##' @param nsplit Number of CV splits. Defaults to 5.
##' @param ... default arguments to covarem().
##' @return List containing (1) the set of coefficients
##' @export
blockcv <- function(cl, folds, cv.gridsize,
                    ylist, countslist, X,
                    mean_lambdas,
                    pie_lambdas,
                    destin,
                    ...){

  ## iimat looks like this:
  ## (#, ialpha, ibeta, ifold)
  ## 1, 1, 1, 1
  ## 2, 1, 1, 2
  ## 3, 1, 1, 3
  ## ...
  ## 1125, 15, 15, 5
  nfold = length(folds)
  iimat = make_iimat(cv.gridsize, nfold)
  iimax = nrow(iimat)
  parallel::parLapplyLB(cl,
  ## lapply(
      1:iimax, function(ii){
    ialpha = iimat[ii,"ialpha"]
    ibeta = iimat[ii,"ibeta"]
    ifold = iimat[ii,"ifold"]
    print('iii=')
    print(c(ialpha, ibeta, ifold))
    one_job(ialpha, ibeta, ifold, folds, destin,
            mean_lambdas, pie_lambdas,
            ylist, countslist, X, ...)
  })
}

##' Synopsis: BLOCKED cross-validation function goes here.

##' CV wrapper for covarem().
##' @param nsplit Number of CV splits. Defaults to 5.
##' @param ... default arguments to covarem().
##' @return List containing (1) the set of coefficients
##' @export
blockcv_fitmodel <- function(cl, folds, destin,
                             ylist, countslist, X,
                             mean_lambdas,
                             pie_lambdas,
                             ...){

  iimat = make_iimat_small(cv.gridsize)
  iimax = nrow(iimat)
  parallel::parLapplyLB(cl,
  ## 1:iimax, function(ii){
  ## lapply(
      1:iimax, function(ii){
    ialpha = iimat[ii, "ialpha"]
    ibeta = iimat[ii, "ibeta"]
    print('iii=')
    print(c(ialpha, ibeta))
    one_job_refit(ialpha, ibeta, folds, destin,
                  mean_lambdas, pie_lambdas,
                  ylist, countslist, X,
                  ...)
  })

}




##' Define the folds.
blockcv_make_folds <- function(ylist, nfold, verbose=FALSE){

  if(verbose) print("Large consecutive time blocks used for CV (block type 1)")
  TT = length(ylist)
  endpoints = round(seq(from = 1, to = TT, length = nfold+1))
  inds = Map(function(a,b){(a+1):b},
             endpoints[-length(endpoints)],
             endpoints[-1])
  return(inds)
}

##' Define the folds for 1-hr-split CVs.
blockcv_hourlong_make_folds <- function(ylist, nfold, verbose=FALSE){

  if(verbose) print("Hour-long time blocks used for CV (block type 2)")
  ## Make hour-long index list
  TT = length(ylist)
  endpoints = round(seq(from = 1, to = TT, by = 20))
  inds = Map(function(a,b){(a+1):b},
             endpoints[-length(endpoints)],
             endpoints[-1])

  ## Further make these into five blocks of test indices.
  test.ii.list = lapply(1:5, function(ifold){
    test.ii = seq(from=ifold, to=length(inds), by=5)
  })
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

  ## Run algorithm on training data,
  res.train = covarem(ylist = train.dat,
                      countslist = train.count,
                      X = train.X,
                      ## refit = FALSE, ## Is this necessary? Think about it for a sec.
                      mean_lambda = mean_lambdas[ibeta],
                      pie_lambda = pie_lambdas[ialpha],
                      ...)
  ## assert_that(!refit)

  ## Assign mn and pie
  pred = predict.covarem(res.train, newx = test.X)
  stopifnot(all(pred$newpie >= 0))

  ## Evaluate on test data, by calculating objective (penalized likelihood)
  cvscore = objective_overall_cov(mu = pred$newmn,
                                  pie = pred$newpie,
                                  sigma = pred$sigma,
                                  ylist = test.dat,
                                  countslist = NULL,
                                  pie_lambda = 0,
                                  mean_lambda = 0,
                                  alpha = res.train$alpha,
                                  beta = res.train$beta)

  ## Save the CV results
  filename = paste0(ialpha, "-", ibeta, "-", ifold, "-cvscore.Rdata")
  save(cvscore, file=file.path(destin, filename))
}



#############################
### Helper to run one job ###
#############################

##' Refit one job
one_job_refit <- function(ialpha, ibeta, folds, destin,
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
##' @param cv.gridsize CV grid size.
##' @param nfold Number of of CV folds.
make_iimat <- function(cv.gridsize, nfold){
  iimat = expand.grid(ialpha = 1:cv.gridsize,
                      ibeta = 1:cv.gridsize,
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
##' @param cv.gridsize CV grid size.
make_iimat_small <- function(cv.gridsize){
  iimat = expand.grid(ialpha = 1:cv.gridsize,
                      ibeta = 1:cv.gridsize)
  iimat = cbind(as.numeric(rownames(iimat)), iimat)
  return(iimat)
}



##' Aggregate results from blocked CV.
##'
##' @param destin destination folder for output.
##' @param cv.gridsize Grid size for cross validation.
##' @param nfold Number of folds.
##'
##' @export
blockcv_aggregate <- function(destin, cv.gridsize, nfold){
  cvscore.array = array(0, dim=c(cv.gridsize, cv.gridsize, nfold))
  cvscore.mat = matrix(0, nrow = cv.gridsize, ncol = cv.gridsize)
  for(ialpha in 1:cv.gridsize){
    for(ibeta in 1:cv.gridsize){
      cvscores <- sapply(1:nfold, function(igrid){
        filename = paste0(ialpha, "-", ibeta, "-", igrid, "-cvscore.Rdata")
        load(file.path(destin, filename))
        return(cvscore)
      })
      cvscore.mat[ialpha, ibeta] = mean(cvscore)
      cvscore.array[ialpha, ibeta, ] = cvscore
    }
  }
  return(list(cvscore.array = cvscore.array,
              cvscore.mat = cvscore.mat))
}
