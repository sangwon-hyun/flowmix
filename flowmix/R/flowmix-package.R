##' flowmix: computing sparse mixture of regressions for flow cytometry data.
##'
##' The main function for conducting analyses is \code{flowmix()}.
##'
##' Try the built-in vignette (`vignette("flowmix")`), for a step-by-step
##' tutorial.
##'
##' @keywords internal 
"_PACKAGE"
##' @name flowmix
##' @useDynLib flowmix, .registration = TRUE
##' @importFrom glmnet glmnet
##' @importFrom glmnet coef.glmnet
##' @import dplyr
##' @importFrom ellipse ellipse
##' @importFrom MASS mvrnorm
##' @importFrom abind abind
##' @useDynLib flowmix
##' @import Rcpp
##' @importFrom Rcpp sourceCpp
##' @importFrom grDevices rgb
##' @import RcppArmadillo
##' @import RcppEigen

# TODO: prevent warning about replacing RcppArmadillo imports by RcppEigen 
# imports

NULL
#> NULL




## This is to fix check()'s global definition errors
pie_lambdas = NULL
cvscore = NULL
objectives = NULL
mean_lambdas = NULL
res = NULL
nfold = NULL
