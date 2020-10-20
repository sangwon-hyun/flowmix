##' Default printing functionality for flowmix() results.
##'
##' @param obj An object of class "flowmix"
##'
##' @return Nothing.
##'
##' @export
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



##' "Prettifies" results.
prettify <- function(res, signif_digit=2){

  ## Reorder clusters in decreasing order of diam
  res <- reorder_clust(res)

  ## Prettify betas
  betamat = do.call(cbind, res$beta)
  all.zero.rows = which(apply(betamat, 1, function(myrow)all(abs(myrow)<1E-8)))
  if(length(all.zero.rows) > 0){
    betamat = betamat[-all.zero.rows,, drop=FALSE]
  }
  betamat = betamat %>% Matrix::Matrix(sparse=TRUE) %>% signif(signif_digit)
  colnames(betamat) = unlist(Map(function(a,b) paste0(a, ", ", b),
                                 paste0("clust-", rep(1:res$numclust, each=res$dimdat)),
                                 colnames(betamat)))

  ## Prettify alphas
  alphamat = res$alpha %>% t %>% Matrix::Matrix(sparse=TRUE) %>% signif(signif_digit)
  all.zero.rows = which(apply(alphamat, 1, function(myrow)all(abs(myrow)<1E-8)))
  if(length(all.zero.rows) > 0){
    alphamat = alphamat[-all.zero.rows,, drop=FALSE]
  }

  ## Return
  return(list(alphamat = alphamat,
              betamat = betamat))
}
