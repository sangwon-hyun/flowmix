## ## Synopsis: contains the parallelized main CV wrapper, and all helper functions
## ## related to cross-validations.

## ##' CV wrapper for covarem().
## ##' @param nsplit Number of CV splits. Defaults to 5.
## ##' @param ... default arguments to covarem().
## ##' @return List containing (1) the set of coefficients
## parallel_cv_test.covarem <- function(ylist, X,
##                                      mean_lambda,
##                                      prob_lambda,
##                                      nsplit = 5,
##                                      numfork = 3,
##                                      verbose = FALSE,
##                                      destin = "~",
##                                      multicore.cv = FALSE,
##                                      warm_start=TRUE,
##                                      ...){
##   ## Create CV split indices
##   assert_that(nsplit >= 2)
##   mysplits = cvsplit(ylist, nsplit = nsplit)

##   ## Save data once
##   save(ylist, X,
##        ## beta_lambdas, alpha_lambdas,
##        ## destin, gridsize, splits, nsplit,
##        file=file.path(destin, "data.Rdata"))

##   ## Parallelize for one pair of lambda values.
##   do_one_pair = function(ylist, X, splits, nsplit, refit, mean_lambda,
##                          prob_lambda, multicore.cv, gridsize, destin, ...){

##       ## The rest is similar to move_to_up() or move_to_left().
##       cvres = get_cv_score(ylist, X, splits, nsplit, refit,
##                            ## Additional arguments for covarem
##                            mean_lambda = mean_lambda
##                            prob_lambda = prob_lambda,
##                            multicore.cv = FALSE,
##                            ...)

##     ## Save a file of the format 1-2-2 where 1 is the number of restart, 2 is
##     ## the fold number. The file will contain two lines (1) the start and end
##     ## time, and (2) the duration (subtract start and end), and (3) Number of
##     ## the iterations it took.

##       ## Get the fitted results on the entire data
##       res = covarem(ylist = ylist, X = X,
##                     mean_lambda = beta_lambdas[ibeta],
##                     prob_lambda = alpha_lambdas[ialpha],
##                     ...)
##       saveres(res = res,
##               cvres = cvres,
##               ialpha = ialpha, ibeta = ibeta, destin = destin,
##               beta_lambdas = beta_lambdas,
##               alpha_lambdas = alpha_lambdas)
##     return(NULL)
##   }

##   do_one_pair(end.ind, ## The rest of the arguments go here
##               ylist, X, mysplits, nsplit, refit, mean_lambdas, prob_lambdas,
##               multicore.cv, gridsize, destin, ...)
## }


##' CV wrapper for covarem().
##' @param nsplit Number of CV splits. Defaults to 5.
##' @param ... default arguments to covarem().
##' @return List containing (1) the set of coefficients
parallel_cv2.covarem <- function(ylist, X,
                                mean_lambdas = NULL,
                                prob_lambdas = NULL,
                                max_mean_lambda = NULL,
                                max_prob_lambda = NULL,
                                gridsize = 9,
                                nsplit = 5,
                                numfork = 3,
                                verbose = FALSE,
                                refit = FALSE,
                                destin = "~",
                                multicore.cv = FALSE,
                                cl=NULL,
                                ...){

  ## Printing some information about the parallelism
  if(verbose==TRUE){
    cat("At most ", numfork * nsplit, " cores will be used.", fill = TRUE)
    cat("Parallel CV output saved to ", destin, fill = TRUE)
  }

  ## Basic checks
  stopifnot(length(mean_lambdas) == length(prob_lambdas))
  assert_that(!is.null(max_mean_lambda) | !is.null(mean_lambdas) )
  assert_that(!is.null(max_prob_lambda) | !is.null(prob_lambdas) )
  if(is.null(mean_lambdas)){
    mean_lambdas = c(exp(seq(from = -8, to = log(max_mean_lambda), length = gridsize)))
  }
  if(is.null(prob_lambdas)){
    prob_lambdas = c(exp(seq(from = -8, to = log(max_prob_lambda), length = gridsize)))
  }

  ## Create CV split indices
  assert_that(nsplit >= 2)
  mysplits = cvsplit(ylist, nsplit = nsplit) ## Too big for my liking; it
                                             ## because hundreds of megabytes
                                             ## with T=4000; but OK sure for
                                             ## now. TODO: make MUCH smaller

  ## If warm starts are not needed, do the following:
    assert_that(!multicore.cv)
    assert_that(!is.null(cl),
                msg=paste0("When you disable warm starts, you must provide a |cl| object!",
                "Making a simple 1-core local cluster is easy: cl=makePSOCKcluster(1)"))

    ## Save data once
    save(ylist, X,
         file=file.path(destin, "data.Rdata"))

    ## Parallelize for one pair of lambda values.
    do_one_pair = function(ind, end.ind,
                           ## The rest of the arguments go here
                           ylist, X, splits, nsplit, refit,
                           beta_lambdas, alpha_lambdas, multicore.cv,
                           gridsize, destin,
                           ...){

      ## Temporary
      start.time = Sys.time()

      ## Redefine which lambda indices correspond to ind in 1:gridsize^2
      ialpha =  ceiling(ind/ gridsize)
      ibeta = (ind-1) %% gridsize + 1

      ## Set seed.
      set.seed(ind)

      ## Check whether this version has been done already.
      if(verbose) cat("ialpha, ibeta are:", ialpha, ibeta, "are attempted.", fill=TRUE)
      already_done = checkres(ialpha, ibeta, destin)
      if(already_done) return(NULL)

      tryCatch({

      ## Get the fitted results on the entire data
      res = covarem(ylist = ylist, X = X,
                    mean_lambda = beta_lambdas[ibeta],
                    prob_lambda = alpha_lambdas[ialpha],
                    ...)

      ## Temporary
      end.time = Sys.time()
      lapsetime = round(difftime(Sys.time(), start.time,
                                 units = "secs"), 0)
      ## End of temporary
      saveres(res = res,
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

    ## Actually do the "brute force" parallelization
    if(verbose){
      cat("Brute force parallelizing on ", length(cl), "cores.", fill=TRUE)
    }

    end.ind = gridsize^2
    parallel::parLapplyLB(cl, end.ind:1, do_one_pair, end.ind,
                ## The rest of the arguments go here
                ylist, X, mysplits, nsplit, refit, mean_lambdas,
                prob_lambdas, multicore.cv, gridsize, destin, ...)

}
