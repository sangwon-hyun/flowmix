## Synopsis: contains the main CV wrapper, and all helper functions related to
## cross-validations.

##' CV wrapper for flowmix().
##'
##' @param nsplit Number of CV splits. Defaults to 5.
##' @param ... default arguments to flowmix().
##'
##' @return List containing (1) the set of coefficients
##'
cv.flowmix <- function(ylist, X, mean_lambdas = NULL,
                       prob_lambdas = NULL,
                       max_mean_lambda = NULL,
                       max_prob_lambda = NULL,
                       gridsize = 5,
                       nsplit = 5,
                       verbose = FALSE,
                       refit = FALSE,
                       mc.cores = 1,
                       ...){

}
