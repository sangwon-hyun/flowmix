##' Synopsis: BLOCKED cross-validation function goes here.


##' CV wrapper for covarem().
##' @param nsplit Number of CV splits. Defaults to 5.
##' @param folds CV folds.
##' @param cv_gridsize Grid size of CV.
##' @param destin Destination directory.
##' @param mean_lambdas List of regularization parameters for mean model.
##' @param pie_lambdas List of regularization parameters for pie model.
##' @param iirange Which of the \code{cv_gridsize^2*nfold*nrep} jobs (e.g. of
##'   the form 5-3-2-5) to run. Defaults to \code{NULL}, in which case all jobs
##'   are run in order of these indices i.e. 1-1-1-1, then 1-1-1-2, and so on.
##' @param ... default arguments to covarem().
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


##' CV wrapper for covarem().
##'
##' @param nsplit Number of CV splits. Defaults to 5.
##' @param ... Other arguments to covarem().
##'
##' @return List containing (1) the set of coefficients
##'
##' @export
blockcv_fitmodel <- function(cl, destin,
                             cv_gridsize,
                             mean_lambdas,
                             pie_lambdas,
                             ## Arguments to covarem()
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
                  ## Arguments to covarem()
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
##' @param ... Rest of arguments for \code{covarem_once()}.
##'
##' @return Nothing is returned. Instead, a file named "1-1-1-1-cvscore.Rdata"
##'   is saved in \code{destin}. (The indices here are ialpha-ibeta-ifold-irep).
##'
one_job <- function(ialpha, ibeta, ifold, irep, folds, destin,
                    mean_lambdas, pie_lambdas,
                    ## The rest that is needed explicitly for covarem()
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
  filename = paste0(ialpha, "-", ibeta, "-", ifold, "-", irep, "-cvscore.Rdata")
  if(sim){ filename = paste0(isim, "-", filename) } ## temporary
  if(file.exists(file.path(destin, filename))){
    if(sim){ ## temporary
      cat('(isim, ialpha, ibeta, ifold, irep)=(',
          sim,",", ialpha,",", ibeta,",", ifold,",", irep,") are already done.", fill=TRUE)
    } else{
      cat("(ialpha, ibeta, ifold, irep) = (",
          ialpha,",", ibeta,",", ifold,",", irep, ") are already done.", fill=TRUE)
    }
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
    args = args[-which(names(args) %in% "nrep")] ## remove |nrep| prior to feeding to covarem_once().
  }


  tryCatch({

    ## New and better do.call() statement:
    ## res.train = do.call(covarem_once, args) ## Old
    argn <- lapply(names(args), as.name)
    names(argn) <- names(args)
    call <- as.call(c(list(as.name("covarem_once")), argn))
    res.train = eval(call, args)

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
    ## Note to self: if the alpha, beta and lambdas are saved, one can recreate
    ## the fitted solutions easily, by creating the means and pies.

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
                          ## The rest that is needed explicitly for covarem_once()
                          ylist, countslist, X,
                          sim=FALSE, isim=1,
                          ...){

  args = list(...)
  nrep = args$nrep
  for(irep in 1:nrep){

    ## Writing file
    filename = paste0(ialpha, "-", ibeta, "-", irep, "-fit.Rdata")
    if(sim){filename = paste0(isim, "-", filename)} ## Temporary
    if(file.exists(file.path(destin, filename))){
      if(sim){
        cat("Refitting for (isim, ialpha, ibeta, irep) = (",
            isim, ",", ialpha, ",",  ibeta, ",", irep, ") is already done.", fill = TRUE)
      } else {
        cat("Refitting for (ialpha, ibeta, irep) = (",
            ialpha,",", ibeta, ",", irep, ") is already done.", fill = TRUE)
      }
      return(NULL)
    } else {

      ## Get the fitted results on the entire data
      args = list(...)
      args$ylist = ylist
      args$countslist = countslist
      args$X = X
      args$pie_lambda = pie_lambdas[ialpha] ## <-----There was a BUG!! pie_lambda[ialpha]
      args$mean_lambda = mean_lambdas[ibeta]
      if("nrep" %in% names(args)) args = args[-which(names(args) %in% "nrep")] ## remove |nrep| prior to feeding

      ## Call the function.
      argn <- lapply(names(args), as.name)
      names(argn) <- names(args)
      call <- as.call(c(list(as.name("covarem_once")), argn))
      ## set.seed(0)
      ## res2 = eval(call, args)
      ## pdf("~/Desktop/somefile.pdf")
      ## plot(res$objectives)
      ## graphics.off()
      ## Sys.sleep(100)

      ## Save the results
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



##' Aggregate results from blocked CV.
##'
##' @param destin destination folder for output.
##' @param cv_gridsize Grid size for cross validation.
##' @param nfold Number of folds.
##'
##' @export
## blockcv_aggregate <- function(destin, cv_gridsize, nfold){
blockcv_aggregate <- function(destin, cv_gridsize, nfold, nrep,
                              save=FALSE, resfile = "all-cvres.Rdata"){

  ## ## Read the meta data (for |nfold|, |cv_gridsize|, |nrep|)
  ## load(file = file.path(destin, 'meta.Rdata'))

  ## Aggregate the results
  cvscore.array = array(NA, dim=c(cv_gridsize, cv_gridsize, nfold, nrep))
  cvscore.mat = matrix(NA, nrow = cv_gridsize, ncol = cv_gridsize)
  for(ialpha in 1:cv_gridsize){
    for(ibeta in 1:cv_gridsize){
      ## cvscores <- sapply(1:nfold, function(igrid){
      obj = matrix(0, nrow=nfold, ncol=nrep)
      for(ifold in 1:nfold){
        for(irep in 1:nrep){
          filename = paste0(ialpha, "-", ibeta, "-", ifold,
                            "-", irep, "-cvscore.Rdata")
          tryCatch({
            load(file.path(destin, filename))
            ## print(filename)
            cvscore.array[ialpha, ibeta, ifold, irep] = cvscore
            obj[ifold, irep] = objectives[length(objectives)]
          }, error = function(e){})
        }
      }

      ## Pick out the CV scores with the biggest objective value
      cvscores = cvscore.array[ialpha, ibeta,,]
      best.cvscores = unlist(sapply(1:nfold, function(ifold){
        cvscores[which.max(obj[,ifold]), ifold]
      }))
      cvscore.mat[ialpha, ibeta] = mean(best.cvscores, na.rm=TRUE)
    }
  }

  ## Clean a bit
  cvscore.mat[which(is.nan(cvscore.mat), arr.ind=TRUE)] = NA
  rownames(cvscore.mat) = signif(pie_lambdas,3)
  colnames(cvscore.mat) = signif(mean_lambdas,3)


  ## Find the minimum
  mat = cvscore.mat
  min.inds = which(mat == min(mat, na.rm=TRUE), arr.ind=TRUE)

  ## Recent addition
  if(save){
    cat("Saving aggregated results to ", file.path(destin, resfile))
    save(cvscore.array,
         cvscore.mat,
         mean_lambdas,
         pie_lambdas,
         min.inds,
         file = file.path(destin, resfile))
  }

  return(list(cvscore.array = cvscore.array,
              cvscore.mat = cvscore.mat,
              mean_lambdas = mean_lambdas,
              pie_lambdas = pie_lambdas,
              min.inds = min.inds))
}




##' Apply the 1SE rule (results are similar to those of get_optimal_info())
##' @param outputdir Location of the output files. e.g. outputdir="~/output/vanilla"
##' @return NULL.
blockcv_onese <- function(destin,
                          gridsize = 12){
## outputdir="~/repos/flowcy/output/vanilla"

  ## Gather the degrees of freedom
  ## dfmat = aggregateres_df(gridsize, destin)
  ## cvscoremat = aggregateres(gridsize, destin)
  cvscoremat = blockcv_aggregate(destin, cv_gridsize, nfold, nrep)$cvscoremat

  ## Add the line
  ind = ind.min = which(cvscoremat == min(cvscoremat, na.rm=TRUE), arr.ind=TRUE)

  ## Add horizontal lines
  load(file=file.path(destin, paste0(ind[1], "-", ind[2], ".Rdata")))
  se = sd(cvres$all)/sqrt(length(cvres$all))

  ## Apply the 1SE rule
  inds = which((cvscoremat <  cvscoremat[ind[1], ind[2]] + se &
                cvscoremat >  cvscoremat[ind[1], ind[2]] - se &
               dfmat < dfmat[ind[1], ind[2]]),
               arr.ind = TRUE)
  if(nrow(inds) > 0){

    ## Pick the minimum CV error that is most parsimonious.
    all.df = apply(inds, 1, function(myind){dfmat[myind[1],myind[2]]})
    inds = inds[which(all.df==min(all.df)),,drop=FALSE]
    all.cvscores = apply(inds, 1, function(myind){cvscoremat[myind[1],myind[2]]})
    ind.min = inds[which.min(all.cvscores),,drop=FALSE]
    assert_that(nrow(ind.min)==1)

  }

  ## Obtain the best parameters and return
  filename = paste0(ind.min[1], "-", ind.min[2], ".Rdata")
  load(file=file.path(destin, filename))
  lambda_alpha = alpha_lambdas[ind.min[1]]
  lambda_beta = beta_lambdas[ind.min[2]]

  return(list(param = c(lambda_alpha, lambda_beta),
              res = res,
              se = se,
              ind.min = ind.min,
              cvscoremat = cvscoremat,
              dfmat = dfmat))
}



##' Helper to aggregate parallelized CV results and obtain degrees of freedom
##' (DF) estimate, saved in |destin|.
##'
##' @param gridsize Size of CV grid.
##' @param destin Location of saved things.
##'
##' @return Matrix containing estimated degrees of freedom.
blockcv_aggregate_df <- function(gridsize, nrep, destin,
                                 save=FALSE, resfile = "all-cvres-df.Rdata"){

  df.array = obj.array = df.alpha.array = df.beta.array = array(NA, dim=c(gridsize, gridsize, nrep))
  df.mat = df.alpha.mat = df.beta.mat = matrix(NA, ncol=gridsize, nrow=gridsize)
  for(ialpha in 1:gridsize){
    for(ibeta in 1:gridsize){

      ## Objective value
      obj = rep(NA, nrep)
      df = df.alpha = df.beta = rep(NA, nrep)
      for(irep in 1:nrep){

        tryCatch({
          ## Load fitted result
          filename = paste0(ialpha, "-", ibeta, "-", irep, "-fit.Rdata")
          load(file.path(destin, filename))

          ## Calculate DF
          df[irep] = do.call(sum, lapply(res$beta, function(mybeta){
            sum(mybeta[-1,]!=0)})) + sum(res$alpha[,-1]!=0)

          df.alpha[irep] = sum(res$alpha[,-1]!=0)

          df.beta[irep] = do.call(sum, lapply(res$beta, function(mybeta){
            sum(mybeta[-1,]!=0)}))

          ## Also calculate objective function
          objectives = res$objectives
          obj[irep] = objectives[length(objectives)]

        }, error = function(e){})
      }
      df.array[ialpha, ibeta, ] = df
      obj.array[ialpha, ibeta,] = obj
      df.alpha.array[ialpha, ibeta, ] = df.alpha
      df.beta.array[ialpha, ibeta, ] = df.beta

      ## Calculate the df of the best model
      if(!all(is.na(obj))){
        ## df.mat[ialpha, ibeta] = df[which.max(obj, na.rm=TRUE)]
        df.mat[ialpha, ibeta] = df[which(obj == min(obj, na.rm = TRUE))]
      }
    }
  }

  ## Assign to new names
  mat = df.mat
  alpha.array = df.alpha.array
  beta.array = df.beta.array

  ## Recent addition
  if(save){
    cat("Saving aggregated results to ", file.path(destin, resfile))
    save(mat,
         alpha.array,
         beta.array,
         df.array,
         file = file.path(destin, resfile))
  }


  ## return(df.mat)
  return(list(mat = mat,
              alpha.array = alpha.array,
              beta.array = beta.array,
              df.array = df.array,
              obj.array = obj.array))
}


##' Helper to aggregate parallelized CV results and obtain the |res| object, all
##' saved in |destin|.
##'
##' @param gridsize Size of CV grid.
##' @param destin Location of saved things.
##'
##' @return List containing, for every (ialpha, ibeta), the "best" result
##'   objects (best in the sense that it had the best likelihood value out of
##'   the |nrep| replicates.)
blockcv_aggregate_res <- function(gridsize, nrep, destin,
                                  save=FALSE, resfile = "best-res.Rdata"){

  ## df.mat = matrix(NA, ncol=gridsize, nrow=gridsize)
  res.list = list()
  for(ialpha in 1:gridsize){
    for(ibeta in 1:gridsize){

      ## Objective value
      obj = rep(NA, nrep)
      ## df = rep(NA, nrep) #
      res.list.inner = list()
      for(irep in 1:nrep){
        filename = paste0(ialpha, "-", ibeta, "-",
                          irep, "-fit.Rdata")
        tryCatch({
          ## Load fitted result
          load(file.path(destin, filename))

          ## Calculate DF
          ## df[irep] = do.call(sum, lapply(res$beta, function(mybeta){
          ##   sum(mybeta[-1,]!=0)})) + sum(res$alpha[,-1]!=0)
          res.list.inner[[irep]] = res
          ## Also calculate objective function
          obj[irep] = res$objectives[length(res$objectives)]
        }, error = function(e){})
      }

      ## Calculate the df of the best model
      if(!all(is.na(obj))){
        ## df.mat[ialpha, ibeta] = df[which.max(obj)    ]
        res.list[[paste0(ialpha, "-", ibeta)]] = res.list.inner[[which.max(obj)]] ## which.min?
      }
    }
  }

  ## Recent addition
  if(save){
    cat("Saving aggregated results to ", file.path(destin, resfile))
    save(res.list, file=file.path(destin, resfile))
  }

  return(res.list)
}





##' An experimental function to visualize the CV results. It calls several other
##' "aggregate" functions in this file.
blockcv_viz <- function(blocktype = 1, datatype = 75, numclust = 5,
                        cv_gridsize = 7,
                        datadir = "~/Dropbox/research/usc/hpc-output", subfolder="",
                        nrep=5){


  ## Load data
  ## destin = file.path(datadir,
  ##                    paste0("blockcv-", blocktype, "-", datatype, "-", numclust))

  if(is.null(subfolder)) subfolder = ""
  destin = file.path(datadir,
                     paste0("blockcv-", blocktype, "-", datatype, "-", numclust),
                     subfolder)


  ## Get the CV results.
  a = blockcv_aggregate(destin, cv_gridsize = cv_gridsize, nfold = 5, nrep = nrep,
                        save=FALSE, resfile = "all-cvres.Rdata")
  cvscore.mat = a$cvscore.mat
  min.inds = a$min.inds
  ## min.inds = c(3,10)

  ## Also get the #nonzero coefficients
  a = blockcv_aggregate_df(gridsize=cv_gridsize, nrep=nrep, destin=destin)
  dfmat = a$mat

  ## Get the refit covarem results
  a = blockcv_aggregate_res(gridsize=cv_gridsize, nrep = nrep, destin=destin)
  res = a[[paste0(min.inds[1] , "-", min.inds[2])]]
  if(is.null(res)) stop(paste0("The model with lambda indices (", min.inds[1], ",", min.inds[2], ") is not available."))

  ## Get the covariates
  X = res$X

  ## Print coefficients
  betalist =  lapply(1:numclust, function(iclust){
    rownames(res$beta[[iclust]])[-1] = colnames(X)
    cf = (res$beta[[iclust]][-1,, drop=FALSE])
    ## cf = round(cf[order(abs(res$beta[[iclust]][-1,]), decreasing=TRUE),,drop=FALSE],3)

    ## Remove the rows that are all zero
    all.zero.rows = which(apply(cf, 1, function(myrow)all(is.na(myrow))))
    if(length(all.zero.rows)>0){
      cf = cf[-all.zero.rows,, drop=FALSE]
    }
    round(Matrix::Matrix(cf, sparse=TRUE),3)
  })
  names(betalist) = paste0("Beta matrix, cluster ", 1:numclust)
  print("Betas are:")
  pretty.betas = betalist
  print(pretty.betas)
  print("Alphas are:")
  colnames(res$alpha)[-1]= colnames(X)
  alpha = t(res$alpha)
  alpha[which(abs(alpha)<1E-5)] = 0
  pretty.alphas = round(Matrix::Matrix(alpha,sparse=TRUE),3)
  print(pretty.alphas)


  ## Get the sigmas
  one_dim = (res$dimdat == 1)
  if(one_dim){
    sigmas = sqrt(res$sigma[,1,])
    names(sigmas) = paste0("Cluster ", 1:numclust)
    pretty.sigmas = sqrt(sigmas)
  } else {
    sigmas = lapply(1:numclust, function(iclust){
      diag(res$sigma[iclust,,])
    })
    names(sigmas) = paste0("Cluster ", 1:numclust)
    pretty.sigmas = lapply(sigmas, sqrt)
  }
  print("Sigmas are:")
  print(pretty.sigmas)

  ## If data only has one dimension, plot it.
  if(one_dim){

    par(mfrow=c(3,1))
    ## Plotting the fitted means and pies:
    load(file.path(destin,"meta.Rdata"))
    TT = length(ylist)
    cols = RColorBrewer::brewer.pal(numclust, "Set2")
    matplot(res$mn[,1,], type='l', lty=1, lwd=1, ylim=range(unlist(ylist)),
            ylab="", xlab="Time", axes=FALSE) ## Base plot
    axis(1);    axis(2)
    title(main="Fitted means and data", cex.main=2)

    ## ## Add demarcation for block boundaries
    ## TT = 308-12
    ## nfold = 5
    ## folds = blockcv_make_folds(1:TT, nfold, verbose = TRUE)
    ## for(iclust in 1:numclust){
    ##   if(iclust%%2==0) abline(v=range(folds[[iclust]]), lwd=3, col="blue")
    ## }

    ## Plotting the fitted means and pies:
    mx = max(unlist(countslist))
    for(tt in 1:TT){
      y = ylist[[tt]]
      ct = countslist[[tt]]/mx
      points(x=rep(tt, length(y)), y=y, col=rgb(0,0,0,ct), pch=15, cex=1)
    }

    ## Replot the means
    for(iclust in 1:numclust){
      lines(res$mn[,1,iclust], type='o', lty=1, lwd=.1, cex=res$pie[,iclust]*5, pch=15, col=cols[iclust])
      lines(res$mn[,1,iclust], type='l', lty=1, lwd=0.1, col=1)
    }
    sigmas = sqrt(res$sigma[,1,])
    sigma_mat = matrix(sqrt(res$sigma[,1,]), ncol=numclust, nrow=TT, byrow=TRUE)
    matlines(res$mn[,1,] + 2*sigma_mat, lty=1, lwd=2, col=cols)
    matlines(res$mn[,1,] - 2*sigma_mat, lty=1, lwd=2, col=cols)

    for(iclust in 1:numclust){
      up = res$mn[,1,iclust] + 2*sigma_mat[,iclust]
      dn = res$mn[,1,iclust] - 2*sigma_mat[,iclust]
      polygon(c(1:length(up),rev(1:length(up))),
              c(up,rev(dn)),col=grDevices::adjustcolor( cols[iclust], alpha.f = 0.2),
              border=NA)
    }
    text(x=rep(-4, numclust), y=res$mn[1,1,], label=paste0("Clust ", 1:numclust), cex=1.5)

    ## Plot the pies as well.
    matplot(res$pie, type='l', lty=1, lwd=3,
            col=cols, ylab="", xlab="Time", ylim=c(0,1), axes=FALSE)
    axis(1);    axis(2)
    abline(h=seq(from=0,to=1, by=0.2), col='grey90', lwd=2, lty=3)
    matlines(res$pie, type='l', lty=1, lwd=3,
            col=cols, axes=FALSE)
    axis(1); axis(2)
    title(main="Cluster probabilities", cex.main=2)

    ## Plot the data
    cols = RColorBrewer::brewer.pal(ncol(X), "Set3")
    matplot(X, type='l', lty=1, ylab="", xlab="Time", lwd=1, col=cols, axes=FALSE)
    axis(1)
    axis(2)
    title(main="Covariates", cex.main=2)
    ## colnames(X)[1:3] = c("b1", "b2", "b2")## temporary
    if(!is.null(colnames(X)))legend("bottomright", lwd=1, col=cols, legend=colnames(X), cex=1.3,
                                    bg="white",  ncol=2)

  }

  return(list(bestres = res,
              cvscore.mat=cvscore.mat,
              min.inds=min.inds,
              dfmatg = dfmat,
              pretty.alphas = pretty.alphas,
              pretty.betas = pretty.betas,
              pretty.sigmas = pretty.sigmas))
              ## mean_lambdas = mean_lambdas,
              ## pie_lambdas = pie_lambdas))
}




##' Temporary function to prettify covarem object res.
prettify <- function(res){

  ## Get the covariates
  X = res$X
  numclust=res$numclust

  ## Print coefficients
  betalist =  lapply(1:numclust, function(iclust){
    rownames(res$beta[[iclust]])[-1] = colnames(X)
    cf = (res$beta[[iclust]][-1,, drop=FALSE])
    ## cf = round(cf[order(abs(res$beta[[iclust]][-1,]), decreasing=TRUE),,drop=FALSE],3)

    ## Remove the rows that are all zero
    all.zero.rows = which(apply(cf, 1, function(myrow)all(myrow==0)))
    if(length(all.zero.rows)>0){
      cf = cf[-all.zero.rows,, drop=FALSE]
    }
    round(Matrix::Matrix(cf, sparse=TRUE),3)
  })
  names(betalist) = paste0("Beta matrix, cluster ", 1:numclust)
  print("Betas are:")
  pretty.betas = betalist
  print(pretty.betas)
  print("Alphas are:")
  colnames(res$alpha)[-1]= colnames(X)
  alpha = t(res$alpha)
  alpha[which(abs(alpha)<1E-5)] = 0
  pretty.alphas = round(Matrix::Matrix(alpha,sparse=TRUE),3)
  print(pretty.alphas)
}



##' Aggregation wrapper, for simulations
blockcv_sim <- function(cvnum = 0, blocktype = 1, datatype = 8, numclust = 2,
                       cv_gridsize = 7,
                       datadir = "~/Dropbox/research/usc/hpc-output", subfolder=""){

  ## Aggregate from CV.

  ## CV indices.

  ## CV indices.

}
