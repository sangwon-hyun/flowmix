## Synopsis: contains the parallelized main CV wrapper, and all helper functions
## related to cross-validations.

##' CV wrapper for covarem().
##' @param nsplit Number of CV splits. Defaults to 5.
##' @param ... default arguments to covarem().
##' @return List containing (1) the set of coefficients
parallel_cv_prebinned.covarem <- function(ylist, X,
                                mean_lambdas = NULL,
                                pie_lambdas = NULL,
                                max_mean_lambda = NULL,
                                max_pie_lambda = NULL,
                                gridsize = 9,
                                nsplit = 5,
                                splits = NULL,
                                numfork = 3,
                                verbose = FALSE,
                                refit = FALSE,
                                destin = "~",
                                cl=NULL,
                                tester = FALSE,
                                train_ybin_list_by_split,
                                train_counts_list_by_split,
                                ylist_orig,
                                ...){

  ## Printing some information about the parallelism
  if(verbose==TRUE){
    cat("Parallel CV output saved to ", destin, fill = TRUE)
    cat("Using ", length(cl), "cores.", fill=TRUE)
  }

  ## Basic checks
  stopifnot(length(mean_lambdas) == length(pie_lambdas))
  assert_that(!is.null(max_mean_lambda) | !is.null(mean_lambdas) )
  assert_that(!is.null(max_pie_lambda) | !is.null(pie_lambdas) )
  if(is.null(mean_lambdas)){
    mean_lambdas = c(exp(seq(from = -8, to = log(max_mean_lambda), length = gridsize)))
  }
  if(is.null(pie_lambdas)){
    pie_lambdas = c(exp(seq(from = -8, to = log(max_pie_lambda), length = gridsize)))
  }
  assert_that(!is.null(splits))
  assert_that(!is.null(cl),
                msg=paste0("You must provide a |cl| object!",
                "Making a simple 1-core local cluster is easy: cl=makePSOCKcluster(1)"))

  ## Parallelize for one pair of lambda values.
  do_one_pair = function(ind, end.ind,
                         ## The rest of the arguments go here
                         ylist, X, splits, nsplit, refit,
                         beta_lambdas, alpha_lambdas,
                         gridsize, destin,
                         train_ybin_list_by_split,
                         train_counts_list_by_split,
                         ylist_orig,
                         ...){

    ## Redefine which lambda indices correspond to ind in 1:gridsize^2
    ialpha =  ceiling(ind/ gridsize)
    ibeta = (ind-1) %% gridsize + 1

    ## Check whether this version has been done already.
    if(verbose) cat("ialpha, ibeta are:", ialpha, ibeta, "are attempted.", fill=TRUE)
    already_done = checkres(ialpha, ibeta, destin)
    if(already_done) return(NULL)

    tryCatch({
    ## The rest is similar to move_to_up() or move_to_left().
    cvres = get_cv_score_prebinned(ylist, X, splits, nsplit, refit,
                                   counts_list,
                                   train_ybin_list_by_split,
                                   train_counts_list_by_split,
                                   ylist_orig,
                                   ## Additional arguments for covarem
                                   mean_lambda = beta_lambdas[ibeta],
                                   pie_lambda = alpha_lambdas[ialpha],
                                   ...)

    ## Get the fitted results on the entire data
    res = covarem(ylist = ylist, X = X,
                  ## Additional arguments for covarem
                  mean_lambda = beta_lambdas[ibeta],
                  pie_lambda = alpha_lambdas[ialpha],
                  ...)

    saveres(res = res,
            cvres = cvres,
            ialpha = ialpha, ibeta = ibeta, destin = destin,
            beta_lambdas = beta_lambdas,
            alpha_lambdas = alpha_lambdas,
            lapsetime = lapsetime)
    }, error = function(err) {
      err$message = paste(err$message,
                          "\n(No file will be saved for the lambdas ",
                          alpha_lambdas[ialpha],", ", beta_lambdas[ibeta],
                          " whose indices are", ialpha,", ", ibeta, " .)",sep="")
      cat(err$message, fill=TRUE)
      warning(err)})
    return(NULL)

  }

  ## temporary feature
  end.ind = gridsize^2
  if(tester){
    lapply(end.ind:1, do_one_pair, end.ind,
           ## The rest of the arguments go here
           ylist, X, splits, nsplit, refit, mean_lambdas, pie_lambdas,
           gridsize, destin,
           train_ybin_list_by_split,
           train_counts_list_by_split,
           ylist_orig,
           ...)
  } else {
    parallel::parLapplyLB(cl, end.ind:1, do_one_pair, end.ind,
                          ## The rest of the arguments go here
                          ylist, X, splits, nsplit, refit, mean_lambdas,
                          pie_lambdas, gridsize, destin,
                          train_ybin_list_by_split,
                          train_counts_list_by_split,
                          ylist_orig,
                          ...)
  }
}

##' Helper to load parallelized CV results, saved in |destin|.
loadres <- function(ialpha, ibeta, destin){
  filename = paste0(ialpha, "-", ibeta, ".Rdata")
  load(file=file.path(destin, filename))
  return(res)
}

##' Helper to save parallelized CV results, saved in |destin|.
saveres <- function(res, cvres, ialpha, ibeta, destin, alpha_lambdas, beta_lambdas,
                    lapsetime## Temporary.
                    ){
  filename = paste0(ialpha, "-", ibeta, ".Rdata")
  save(res, cvres, alpha_lambdas, beta_lambdas,
       lapsetime, ## Temporary.
       file=file.path(destin, filename))
}

##' Helper to see if CV results for (ialpha, ibeta), saved in |destin|, have
##' already been run (i.e. by seeing if the file exists)
checkres <- function(ialpha, ibeta, destin){
  filename = paste0(ialpha, "-", ibeta, ".Rdata")
  return(file.exists(file=file.path(destin, filename)))
}


##' Helper to aggregate parallelized CV results, saved in |destin|.
##' @param gridsize Size of CV grid.
##' @param destin Location of saved things.
##' @param numfold Number of CV folds.
##' @return gridsize x gridsize Matrix containing average CV scores.
aggregateres <- function(gridsize, destin){

  cvscoremat = matrix(NA, gridsize, gridsize)
  for(ialpha in 1:gridsize){
    for(ibeta in 1:gridsize){
      filename = paste0(ialpha, "-", ibeta, ".Rdata")

      ## Check that results exist
      file_exists = file.exists(file = file.path(destin, filename))
      ## assert_that(file_exists)
      if(!file_exists){
        cat(filename, "doesn't exist.", fill=TRUE)
      } else {
        ## Load CV score and insert in matrix
        load(file = file.path(destin, filename))
        cvscoremat[ialpha, ibeta] = cvres$mean
      }
    }
  }
  return(cvscoremat)
}

##' Helper to aggregate parallelized CV results and obtain degrees of freedom
##' (DF) estimate, saved in |destin|.
##' @param gridsize Size of CV grid.
##' @param destin Location of saved things.
##' @return Matrix containing estimated degrees of freedom.
aggregateres_df <- function(gridsize, destin){

  dfmat = matrix(NA, gridsize, gridsize)
  for(ialpha in 1:gridsize){
    for(ibeta in 1:gridsize){
      filename = paste0(ialpha, "-", ibeta, ".Rdata")

      file_exists = file.exists(file = file.path(destin, filename))
      if(!file_exists){
        cat(filename, "doesn't exist.", fill=TRUE)
      } else {
      ## assert_that(file.exists(file = file.path(destin, filename)))
      load(file = file.path(destin, filename))
      mydf = do.call(sum, lapply(res$beta, function(mybeta){
        sum(mybeta[-1,]!=0)})) + sum(res$alpha[,-1]!=0)
      dfmat[ialpha, ibeta] = mydf
      }
    }
  }
  return(dfmat)
}


##' Same as \code{aggregateres()}, but taking all folds!
##' @param gridsize Size of CV grid.
##' @param destin Location of saved things.
##' @param numfold Number of CV folds.
##' @return Array containing /all/ CV scores. Dimension is gridsize x gridsize x
##'   numfold.
aggregateres_allfolds <- function(gridsize, destin, numfold){

  cvscorearray = array(NA, c(gridsize, gridsize, numfold))
  for(ialpha in 1:gridsize){
    for(ibeta in 1:gridsize){
      filename = paste0(ialpha, "-", ibeta, ".Rdata")
      ## Check that results exist
      ## cat("(", ialpha, ibeta, ")", fill=TRUE)
      assert_that(file.exists(file = file.path(destin, filename)))
      ## Load CV score and insert in matrix
      load(file = file.path(destin, filename))
      cvscorearray[ialpha, ibeta, ] = cvres$all
    }
  }
  return(cvscorearray)
}

##' Gets all the lambdas, from one of the results.
##' @param destin Directory containing output.
##' @return List containing the regularization parameter values used for CV.
getlambdasres <- function(destin){
  ialpha = ibeta = 1
  filename = paste0(ialpha, "-", ibeta, ".Rdata")
  assert_that(file.exists(file = file.path(destin, filename)))
  ## Load CV score and insert in matrix
  load(file = file.path(destin, filename))
  return(list(alpha = alpha_lambdas,
              beta = beta_lambdas))
}

##' Gets the fitted object at (ialpha, ibeta).
##' @param ialpha Number between 1-gridsize.
##' @param ibeta Number between 1-gridsize.
##' @param destin Directory containing output.
getres <- function(ialpha, ibeta, destin){
  filename = paste0(ialpha, "-", ibeta, ".Rdata")
  assert_that(file.exists(file = file.path(destin, filename)))
  ## Load CV score and insert in matrix
  load(file = file.path(destin, filename))
  return(res)
}


##' Helper to get optimal model's information from (parallelized) CV results.
##' @param gridsize Size of CV grid.
##' @param isim The simulation number.
##' @param outputdir directory containing output of simulation.
##' @param excludecorner temporary feature to exclude the most regularized
##'   result.
##' @return Array containing /all/ CV scores. Dimension is gridsize x gridsize x
##'   numfold.
get_optimal_info <- function(isim=NULL, outputdir="~/repos/flowcy/output/vanilla",
                             gridsize=12, excludecorner=FALSE){
  if(!is.null(isim)){
    destin = file.path(outputdir, paste0("isim-", isim))
  } else {
    destin = outputdir
  }
  ## TODO: eventually, get rid of the isim altogether. The outputdir that
  ## contains all the "3-5.Rdata" type files should be provided.

  ## Collect CV score matrix
  cvscoremat = aggregateres(gridsize, destin)
  if(excludecorner) cvscoremat[gridsize,gridsize] = Inf
  cvscoremat[which(is.na(cvscoremat), arr.ind=TRUE)] = Inf

  ## Get maximizer
  ind = which(cvscoremat == min(cvscoremat), arr.ind=TRUE)
  filename = paste0(ind[1], "-", ind[2], ".Rdata")
  load(file=file.path(destin, filename))
  lambda_alpha = alpha_lambdas[ind[1]]
  lambda_beta = beta_lambdas[ind[2]]

  ## Also calculate the minimum index
  ind.min = which(cvscoremat==min(cvscoremat, na.rm=TRUE), arr.ind=TRUE)

  return(list(param = c(lambda_alpha, lambda_beta),
              res = res,
              ind.min = ind.min,
              alpha_lambdas = alpha_lambdas,
              beta_lambdas = beta_lambdas,
              cvscoremat = cvscoremat))
}





##' Helper function to get the "cl" object, which is an object of class
##' c("SOCKcluster", "cluster").
##' @param numcores NULL, in which case it identifies (almost) /all/ the cores
##'   available to use, and spans them.
##' @return An object of class "cluster" and "SOCKcluster".
get_cl <- function(numcores=NULL){

  ## Obtain all the cores from all the machines available to the master node,
  ## and span them.
  if(is.null(numcores)){
    nodeNames = strsplit(system("scontrol show hostname $SLURM_NODELIST | paste -d, -s", intern=TRUE),
                           ",")[[1]]
    mycommand = "sinfo -o '%40N %c' --Node --long | grep "
    numCores = sapply(nodeNames, function(mynodename){
      mystring = system(paste0(mycommand, mynodename), intern=TRUE)
      return(as.numeric(strsplit(mystring, "\\s+")[[1]][2]))
    })
    machines = do.call(c, unname(Map(function(a,b){rep(a, b-1 )}, nodeNames, unname(numCores))))
    cat("Identifying and allocating", numCores, "cores on",
        length(nodeNames), "nodes (computers) named", nodeNames, fill = TRUE)
    ## Print setting into output stream /and/ a file
    cl = parallel::makePSOCKcluster(machines)

  ## Just span a few cores on the current machine.
  } else {
    assert_that(is.numeric(numcores))
    cl = parallel::makePSOCKcluster(numcores)
  }
  return(cl)
}


##' Making an nsplit length list of binned datasets.
bin_train_splits <- function(ylist, splits, dat.grid, verbose=FALSE, mc.cores = 1){

  ## Setup
  TT = length(ylist)
  nsplit = length(splits[[1]])

  ## For each test/train split, bin the TRAIN data.
  train_ybin_list_by_split = train_counts_list_by_split = list()
  for(test.isplit in 1:nsplit){
    printprogress(test.isplit, nsplit, "splits", fill=TRUE)
    ylist.train = lapply(1:TT, function(tt){
      ind = splits[[tt]][[test.isplit]]
      ylist[[tt]][-ind,]
    })
    obj = bin_many_cytograms(ylist.train, dat.grid, verbose = verbose, mc.cores=mc.cores)
    train_ybin_list_by_split[[test.isplit]] <- obj$ybin_list
    train_counts_list_by_split[[test.isplit]] <- obj$counts_list
  }
  return(list(train_ybin_list_by_split = train_ybin_list_by_split,
              train_counts_list_by_split = train_counts_list_by_split))
}
