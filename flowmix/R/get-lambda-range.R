##' Estimate maximum lambda values by brute force.  First starts with a large
##' initial value \code{max_mean_lambda} and \code{max_prob_lambda}, and runs
##' the EM algorithm on decreasing set of values (sequentially halved). This
##' stops once you see any non-zero coefficients, and returns the *smallest*
##' regularization (lambda) value pair that gives full sparsity. Note that the
##' \code{zero_stabilize=TRUE} option is used in \code{flowmix()}, which
##' basically means the EM algorithm runs only until the zero pattern
##' stabilizes.
##'
##' @param ylist List of responses.
##' @param countslist Multiplicity for particles in \code{ylist}.
##' @param X Covariates.
##' @param numclust Number of clusters.
##' @param max_mean_lambda Defaults to 4000.
##' @param max_prob_lambda Defaults to 1000.
##' @param iimax Maximum value of \eqn{x} for \eqn{2^{-x}} factors to try.
##' @param verbose TRUE for loudness.
##' @param ... Other arguments to \code{flowmix_once()}.
##' @return list containing the two maximum values to use.
##' @examples
##' \dontrun{
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
##' numclust = 4
##' maxres = calc_max_lambda(ybin_list, counts_list, X, numclust, verbose=TRUE,
##'                             nrep = 4,
##'                             ## Function settings
##'                             parallelize = FALSE,
##'                             iimax = 20,
##'                             niter = 1000,
##'                             max_prob_lambda = 10000,
##'                             tol = 1E-3 ## This doesn't need to be so low here.
##'                             )
##'
##' }
##' @export
calc_max_lambda <- function(ylist, countslist = NULL, X, numclust,
                            max_mean_lambda = 4000,
                            max_prob_lambda = 1000,
                            verbose = FALSE,
                            iimax = 16,
                            ...){

  ## Get range of regularization parameters.
  facs = sapply(1:iimax, function(ii) 2^(-ii+1)) ## DECREASING order
  print("running the models once")
  for(ii in 1:iimax){

    ## print_progress(ii, iimax, "regularization values", fill = TRUE)
    cat("###############################################################", fill=TRUE)
    cat("#### lambda_alpha = ", max_prob_lambda * facs[ii],
        " and lambda_beta = ", max_mean_lambda * facs[ii], "being tested.  ", fill=TRUE)
    cat("###############################################################", fill=TRUE)

    res = flowmix_once(ylist = ylist,
                       countslist = countslist,
                       X = X,
                       numclust = numclust,
                       prob_lambda = max_prob_lambda * facs[ii],
                       mean_lambda = max_mean_lambda * facs[ii],
                       verbose = verbose,
                       zero_stabilize = TRUE,
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
        stop(paste0("Max lambdas: ", max_mean_lambda, " and ",
                    max_prob_lambda,
                    " were too small as maximum reg. values. Go up and try again!!"))

      ## If there are *any* nonzero values, return the immediately preceding
      ## lambda values -- these were the smallest values we had found that gives
      ## full sparsity.
      } else {
        ## Check one more time whether the model was actually zero, by fully running it;
        res = flowmix_once(ylist = ylist,
                           countslist = countslist,
                           X = X,
                           numclust = numclust,
                           prob_lambda = max_prob_lambda * facs[ii],
                           mean_lambda = max_mean_lambda * facs[ii],
                           zero_stabilize = FALSE,
                           ...)
        toler = 0
        sum_nonzero_alpha = sum(res$alpha[,-1] > toler)
        sum_nonzero_beta = sum(unlist(lapply(res$beta, function(cf){ sum(cf[-1,] > toler) })))

        ## If there are *any* nonzero values, do one of the following
        if(sum_nonzero_alpha + sum_nonzero_beta != 0){
          return(list(beta = max_mean_lambda * facs[ii-1],
                      alpha = max_prob_lambda *facs[ii-1]))
        }
        ## Otherwise, just proceed to the next iteration.
      }
    }
    cat(fill=TRUE)
  }
}

##' A wrapper for \code{calc_max_lambda}. Saves the two maximum lambda values in
##' a file.
##'
##' @param destin Where to save the output (A two-lengthed list called
##'   "maxres").
##' @param maxres_file Filename for output. Defaults to maxres.Rdata.
##' @param maxdev Maximum value of lambda.
##' @param ... Additional arguments to \code{flowmix()}.
##' @inheritParams calc_max_lambda
##'
##' @return List of the pair of lambda regularization parameter values.
##'
##' @export
get_max_lambda <- function(destin, maxres_file = "maxres.Rdata",
                           ylist,
                           countslist,
                           X,
                           numclust,
                           maxdev,
                           max_prob_lambda,
                           max_mean_lambda,
                           ...){

  if(file.exists(file.path(destin, maxres_file))){
    load(file.path(destin, maxres_file))
    cat("Maximum regularization values are loaded.", fill=TRUE)
    return(maxres)
  } else {
    print(Sys.time())
    cat("Maximum regularization values being calculated.", fill = TRUE)
    cat("with initial lambdas values (alpha and beta):", fill = TRUE)
    print(c(max_prob_lambda, max_mean_lambda));
    maxres = calc_max_lambda(ylist = ylist,
                             countslist = countslist,
                             X = X,
                             numclust = numclust,
                             maxdev = maxdev,
                             ## This function's settings
                             max_prob_lambda = max_prob_lambda,
                             max_mean_lambda = max_mean_lambda,
                             ...)
    save(maxres, file = file.path(destin, maxres_file))
    cat("file was written to ", file.path(destin, maxres_file), fill=TRUE)
    cat("maximum regularization value calculation done.", fill = TRUE)
    print(Sys.time())
    return(maxres)
  }
}



##' Helper function to logarithmically space out R.  \code{length} values linear
##' on the log scale from \code{max} down to \code{min}.
##'
##' @param max Maximum value.
##' @param min Minimum value.
##' @param length Length of the output string.
##' @param min.ratio Factor to multiply to \code{max}.
##'
##' @return Log spaced
##'
##' @export
logspace <- function(max, min=NULL, length, min.ratio = 1E-4){
  if(is.null(min)) min = max * min.ratio
  vec = 10^seq(log10(min), log10(max), length = length)
  stopifnot(abs(vec[length(vec)] - max) < 1E10)
  return(vec)
}
