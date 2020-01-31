##' Estimate maximum lambda values by brute force.  First starts with a large
##' initial value \code{max_lambda_beta} and \code{max_lambda_alpha}, and runs
##' the EM algorithm on decreasing grid of values. Either stop once you see any
##' non-zero coefficients, or run algorithms on the entire grid (Set
##' \code{parallelize=TRUE} and provide a \code{cl} object) and filter this grid
##' (large factor to small), in order to return the *smallest* regularization
##' (lambda) value pair that gives full sparsity. Note that the
##' \code{zero_stabilize=TRUE} option is used in \code{covarem()}, which
##' basically means the EM algorithm runs only until the zero pattern
##' stabilizes.
##'
##' @param ylist List of responses.
##' @param X Covariates.
##' @param numclust Number of clusters.
##' @param max_lambda_beta Defaults to 4000.
##' @param max_lambda_beta Defaults to 1000.
##' @param iimax Maximum value of x for 2^{-x} factors to try.
##' @param parallelize TRUE if jobs are to be run in parallel.
##' @param ... Other arguments to \code{covarem_once()}.
##' @return list containing the two maximum values to use.
##' @examples
##' ## Generate and bin data
##' obj = generate_data_generic(p=5, TT=50, fac=1, nt=7000, dimdat=3)
##' ylist = obj$ylist
##' X = obj$X
##' dat.gridsize = 50
##' dat.grid = make_grid(ylist, gridsize = dat.gridsize)
##' obj = bin_many_cytograms(ylist, dat.grid, mc.cores=4, verbose=TRUE)
##' ybin_list = obj$ybin_list
##' counts_list = obj$counts_list
##'
##' ## No parallel:
##' numclust = 4
##' maxres = calc_max_lambda(ybin_list, counts_list, X, numclust, verbose=TRUE,
##'                             nrep = 4,
##'                             ## Function settings
##'                             parallelize = FALSE,
##'                             iimax = 20,
##'                             niter = 1000,
##'                             max_lambda_alpha = 10000,
##'                             tol = 1E-3 ## This doesn't need to be so low here.
##'                             )
##'
##'
##' # Yes parallel:
##' cl = get_cl(3)
##' parallel::clusterExport(cl, ls())
##' parallel::clusterCall(cl, function(){ load_all("~/repos/flowcy/flowcy")}) ## directory that contains the R package.
##' maxres = calc_max_lambda(ybin_list, counts_list, X, numclust, verbose=TRUE,
##'                             nrep = 4,
##'                             ## Function settings
##'                             parallelize = TRUE,
##'                             cl = cl,
##'                             iimax = 6,
##'                             niter = 1000,
##'                             max_lambda_alpha = 10000,
##'                             tol = 1E-3 ## This doesn't need to be so low here.
##'                             )
calc_max_lambda <- function(ylist, countslist = NULL, X, numclust,
                               max_lambda_beta = 4000,
                               max_lambda_alpha = 1000,
                               verbose=FALSE,
                               iimax = 16,
                               parallelize = FALSE,
                               cl = NULL,
                               ...){

  ## Get range of regularization parameters.
  ## res0 = covarem_getrange(ylist=ylist, X=X, numclust=numclust, niter=2)

  ## ################################
  ## ## First option: mclapply() ####
  ## ################################
  mc.cores = 1
  facs = sapply(1:iimax, function(ii) 2^(-ii+1)) ## DECREASING order
  print("running the models once")
  if(!parallelize){
    for(ii in 1:iimax){
      if(verbose) printprogress(ii, iimax, fill=TRUE)
      cat("============================================================", fill=TRUE)
      cat("lambda_alpha = ", max_lambda_alpha * facs[ii],
          " and lambda_beta = ", max_lambda_beta * facs[ii], "being tested.", fill=TRUE)
      res = covarem_once(ylist = ylist,
                         countslist = countslist,
                         X = X,
                         numclust = numclust,
                         pie_lambda = max_lambda_alpha * facs[ii],
                         mean_lambda = max_lambda_beta * facs[ii],
                         verbose = TRUE,
                         zero_stabilize =TRUE,
                         ...)

      ## Check zero-ness
      toler = 0
      sum_nonzero_alpha = sum(res$alpha[,-1] > toler)
      sum_nonzero_beta = sum(unlist(lapply(res$beta, function(cf){ sum(cf[-1,] > toler) })))


      ## If there are *any* nonzero values, do one of the following
      if(sum_nonzero_alpha + sum_nonzero_beta != 0){


        ## If there are *any* nonzero values at the first iter, prompt a restart
        ## with higher initial lambda values.
        if(ii==1){
          stop(paste0("Max lambdas: ", max_lambda_beta, " and ",
                      max_lambda_alpha,
                      " were too small as maximum reg. values. Go up and try again!!"))


        ## If there are *any* nonzero values, return the immediately preceding
        ## lambda values -- these were the smallest values we had found that gives full sparsity.
        } else {
          return(list(beta = max_lambda_beta * facs[ii-1], alpha = max_lambda_alpha *facs[ii-1]))
        }
      }
    }
  } else {
    assert_that(!is.null(cl))
    reslist = parallel::parLapplyLB(cl, 1:iimax, function(ii){
      if(verbose) printprogress(ii, iimax, fill=TRUE)
      res = covarem(ylist = ylist,
                    countslist = countslist,
                    X = X,
                    numclust = numclust,
                    pie_lambda = max_lambda_alpha * facs[ii],
                    mean_lambda = max_lambda_beta * facs[ii],
                    verbose=TRUE,
                    ...)
      return(res[c("alpha", "beta", "mean_lambda", "pie_lambda")])
    })

    ## Then, filter the grid (large factor to small), to return the smallest
    ## regularization value pair that gives full sparsity.
    print("filtering the results")
    allzero = rep(NA, iimax)
    for(ii in 1:iimax){
      sum_nonzero_alpha = sum(reslist[[ii]]$alpha[,-1] > toler)
      sum_nonzero_beta = sum(unlist(lapply(reslist[[ii]]$beta, function(cf){ sum(cf[-1,] > toler) })))
      allzero[ii] = print(sum_nonzero_alpha + sum_nonzero_beta == 0)
    }
    if(any(allzero)){
      ## Get the SMALLEST coeff such that all zero coefficients were found, and return.
      myfac = facs[max(which(allzero))]
      ## return(c(max_lambda_beta * myfac, max_lambda_alpha * myfac))
      return(list(beta = max_lambda_beta * myfac, alpha = max_lambda_alpha * myfac,
                  reslist = reslist)) ## Addition
    } else {
      stop(paste0("Max lambdas: ", max_lambda_beta, " and ", max_lambda_alpha,
                  " were too small as maximum reg. values. Go up and try again!!"))
    }
  }
}

##' A wrapper for \code{calc_max_lambda}. Saves the two maximum lambda values in
##' a file.
##'
##' @param destin Where to save the output (A two-lengthed list called
##'   "maxres").
##' @param maxres_file Filename for output. Defaults to maxres.Rdata.
##' @param ... Additional arguments to \code{covarem()}.
##'
##' @return No return
get_max_lambda <- function(destin, maxres_file = "maxres.Rdata",
                           ylist,
                           countslist,
                           X,
                           numclust,
                           maxdev,
                           max_lambda_alpha,
                           max_lambda_beta,
                           ...){

  maxres_file = "maxres.Rdata"
  if(file.exists(file.path(destin, maxres_file))){
    load(file.path(destin, maxres_file))
    cat("Maximum regularization values are loaded.", fill=TRUE)
    return(maxres)
  } else {
    print(Sys.time())
    cat("Maximum regularization values being calculated.", fill=TRUE)
    print("with initial lambdas values (alpha and beta):")
    print(c(max_lambda_alpha, max_lambda_beta));
    maxres = calc_max_lambda(ylist = ylist,
                             countslist = countslist,
                             X = X,
                             numclust = numclust,
                             maxdev = maxdev,
                             ## This function's settings
                             parallelize = FALSE,
                             ## cl = cl,
                             max_lambda_alpha = max_lambda_alpha,
                             max_lambda_beta = max_lambda_beta,
                             ...)
    save(maxres, file=file.path(destin, maxres_file))
    cat("maximum regularization value calculation done.", fill=TRUE)
    print(Sys.time())
    return(maxres)
  }
}
