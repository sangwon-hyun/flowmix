##' flowmix: computating sparse mixture of regressions for flow cytometry data.
##'
##' The main function for conducting analyses is \code{flowmix()}.
##'
##' Try the built-in vignette (`vignette("flowmix")`), for a step-by-step
##' tutorial.
##'
##' @docType package
##' @name flowmix
##' @useDynLib flowmix, .registration = TRUE
##' @importFrom glmnet glmnet
##' @importFrom glmnet coef.glmnet
##' @import CVXR
##' @import dplyr
##' @importFrom ellipse ellipse
##' @importFrom sn rsn
##' @importFrom MASS mvrnorm
##' @importFrom abind abind
##' @useDynLib flowmix
##' @import Rcpp
##' @importFrom Rcpp sourceCpp
##' @import RcppArmadillo
##' @import RcppEigen
NULL
#> NULL




## This is to fix check()'s global definition errors
pie_lambdas = NULL
cvscore = NULL
objectives = NULL
mean_lambdas = NULL
res = NULL
nfold = NULL
