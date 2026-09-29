##' Create random weight neural network hidden nodes 
##' 
##' @description Huang, Zhu, and Siew (2006)
##' 
##' @param X An n x p design matrix where p is the 
##' number of covariates (without a column of all 1s 
##' corresponding to an intercept term).
##' @param n.h Number of hidden nodes to output
##' @param a Hidden layer weights are drawn from 
##' Unif(-a, a)
##' 
##' @return An n x n.h matrix of hidden layer outputs.
##' 
##' @export
make_hidden_nodes <- function(X, n.h, a) {
    X_tilde <- cbind(1, X)
    p_tilde <- ncol(X_tilde)
    W <- stats::runif(p_tilde * n.h, -a, a) %>% 
        matrix(nrow = p_tilde)
    
    out_layer <- stats::plogis(X_tilde %*% W)
    return(out_layer)
}