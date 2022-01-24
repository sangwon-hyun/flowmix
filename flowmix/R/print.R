##' Default printing functionality for flowmix() results.
##'
##' @param x An object of class "flowmix".
##' @param ... Remaining parameters.
##'
##' @return Nothing.
##'
##' @export
print.flowmix <- function(x, ...){
  betas = do.call(cbind, lapply(x$beta, function(beta)beta[-1,]))
  alpha =x$alpha[ ,-1]
  print(paste0(sum(betas==0), " out of ", length(betas),
               " regression coefficients for the cluster centers, are zero."))
  print(paste0(sum(alpha==0), " out of ", length(alpha),
               " regression coefficients for the cluster probabilities, are zero."))
  if(x$niter == x$final.iter){
    print(paste0("The algorithm has NOT converged in ", x$niter,
                 "EM iterations."))
  } else {
    print(paste0("The algorithm converged in ", x$final.iter, " EM iterations."))
  }
}
