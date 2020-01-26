##' Estimate maximum lambda values by brute force.  First starts with a large
##' initial value \code{max_lambda_beta} and \code{max_lambda_alpha}, and runs
##' the EM algorithm on decreasing grid of values. Either stop once you see any
##' non-zero coefficients, or run algorithms on the entire grid (Set
##' \code{parallelize=TRUE} and provide a \code{cl} object) and filter this grid
##' (large factor to small), in order to return the *smallest* regularization
##' (lambda) value pair that gives full sparsity.
##'
##' @param ylist List of responses.
##' @param X Covariates.
##' @param numclust Number of clusters.
##' @param max_lambda_beta Defaults to 4000.
##' @param max_lambda_beta Defaults to 1000.
##' @param iimax Maximum value of x for 2^{-x} factors to try.
##' @param parallelize TRUE if jobs are to be run in parallel.
##' @param zero_stabilize Defaults to TRUE, in which case the EM is only run
##'   until the zero pattern in the coefficients stabilize.
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
##' maxres = get_max_lambda(ybin_list, counts_list, X, numclust, verbose=TRUE,
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
##' maxres = get_max_lambda(ybin_list, counts_list, X, numclust, verbose=TRUE,
##'                             nrep = 4,
##'                             ## Function settings
##'                             parallelize = TRUE,
##'                             cl = cl,
##'                             iimax = 6,
##'                             niter = 1000,
##'                             max_lambda_alpha = 10000,
##'                             tol = 1E-3 ## This doesn't need to be so low here.
##'                             )
get_max_lambda <- function(ylist, countslist = NULL, X, numclust,
                               max_lambda_beta = 4000,
                               max_lambda_alpha = 1000,
                               verbose=FALSE,
                               iimax = 16,
                               parallelize = FALSE,
                               cl = NULL,
                               zero_stabilize = TRUE, ## temporary
                               ...){

  ## Get range of regularization parameters.
  ## res0 = covarem_getrange(ylist=ylist, X=X, numclust=numclust, niter=2)

  ## ################################
  ## ## First option: mclapply() ####
  ## ################################
  mc.cores = 1
  facs = sapply(1:iimax, function(ii) 2^(-ii)) ## DECREASING order
  print("running the models once")
  if(!parallelize){
    for(ii in 1:iimax){
      if(verbose) printprogress(ii, iimax, fill=TRUE)
      res = covarem_once(ylist = ylist,
                    countslist = countslist,
                    X = X,
                    numclust = numclust,
                    pie_lambda = max_lambda_alpha * facs[ii],
                    mean_lambda = max_lambda_beta * facs[ii],
                    verbose = TRUE,
                    zero_stabilize = zero_stabilize, ## temporary
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

  ## ## #################################################
  ## ## ## New: Do this on a 2d grid; using parLapply. ##
  ## ## #################################################
  ## res0 = list(max_lambda_beta = 2000, max_lambda_alpha = 20)
  ## facs = sapply(1:8, function(ii) 2^(-ii))

  ## ## Do the full 2d thing.
  ## iimat = expand.grid(1:8, 1:8)
  ## iimax = nrow(iimat)
  ## cl1 = get_cl()
  ## reslist = parLapplyLB(cl1, 1:iimax, function(ii){
  ##   fac1 = 2000 * facs[iimat[ii,1]]
  ##   fac2 = 20 * facs[iimat[ii,2]]
  ##   res = covarem(ylist = ylist,
  ##                 X = X,
  ##                 numclust = numclust,
  ##                 mean_lambda = res0$max_lambda_beta,
  ##                 pie_lambda = res0$max_lambda_alpha,
  ##                 ...)
  ##   return(res[c("alpha", "beta", "mean_lambda", "pie_lambda")])
  ## })

  ## ## Then, filter the grid for the number of zeros.
  ## tol = 1E-8
  ## zeromat = matrix(0, nrow = iimax, ncol = 2)
  ## for(ii in 1:iimax){
  ##   sum_zero_alpha = sum(reslist[[ii]]$alpha <= tol)
  ##   sum_zero_beta = sum(unlist(lapply(res$beta, function(cf){ sum(cf[-1,]==0) })) <= tol)
  ##   zeromat[ii,] = c(sum_zero_alpha, sum_zero_beta)
  ## }

  ## allzero = apply(zeromat, 1, function(myrow){
  ##   all(myrow==0)
  ## })

  ## if(any(allzero)){
  ##   ## Randomly choose a guy out of these, and return.
  ##   inds = iimat[sample(which(allzero), 1),]
  ##   return(c(2000 * facs[inds[1]], 20 * facs[inds[2]]))
  ## }
}
