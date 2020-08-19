##' Default printing functionality for flowmix() results.
##' @param obj An object of class "flowmix"
##' @return Nothing.
print.flowmix <- function(obj){
  betas = do.call(cbind, lapply(obj$beta, function(beta)beta[-1,]))
  alpha =obj$alpha[ ,-1]
  print(paste0(sum(betas==0), " out of ", length(betas),
               " regression coefficients for the cluster centers, are zero."))
  print(paste0(sum(alpha==0), " out of ", length(alpha),
               " regression coefficients for the cluster probabilities, are zero."))
  if(obj$niter == obj$final.iter){
    print(paste0("The algorithm has NOT converged in ", obj$niter,
                 "EM iterations."))
  } else {
    print(paste0("The algorithm converged in ", obj$final.iter, " EM iterations."))
  }
}
