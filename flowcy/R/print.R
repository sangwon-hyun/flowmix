##' Default printing functionality for covarem() results.
##' @param obj An object of class "covarem"
##' @return Nothing.
print.covarem <- function(obj){
  betas = do.call(cbind, lapply(res$beta, function(beta)beta[-1,]))
  alpha =res$alpha[ ,-1]
  print(paste0(sum(betas==0), " out of ", length(betas),
               " regression coefficients for the cluster centers, are zero."))
  print(paste0(sum(obj$alpha==0), " out of ", length(obj$alpha),
               " regression coefficients for the cluster probabilities, are zero."))
  if(obj$niter == obj$final.iter){
    print(paste0("The algorithm has NOT converged in ", obj$niter,
                 "EM iterations."))
  } else {
    print(paste0("The algorithm converged in ", obj$final.iter, " EM iterations."))
  }
}
