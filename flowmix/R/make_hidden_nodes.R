#' Create a (p + 1) x n.h matrix of hidden layer weights 
#' and biases, where all elements are iid Unif(-a, a)
#' @noRd
make_unif_W <- function(p, n.h, a) {
    W <- stats::runif((p + 1) * n.h, -a, a)
    W <- matrix(W, nrow = p + 1)
    return(W)
}

#' Element-wise standard logistic function
#' @noRd 
logistic <- function(V) {
    return(1 / (1 + exp(-V)))
}

#' Create extreme learning machine (ELM) hidden nodes 
#' 
#' @description Huang, Zhu, and Siew (2006)
#' 
#' @param X An n x p design matrix where p is the 
#' number of covariates (X cannot have a column of all 1s 
#' corresponding to an intercept term)
#' @param n.h Number of hidden nodes to output
#' @param a Hidden layer weights are drawn from 
#' Unif(-a, a)
#' 
#' @return An n x n.h matrix resulting from inputting 
#' the data X through the ELM.
#' 
#' @export
make_hidden_nodes <- function(X, n.h, a) {
    p <- ncol(X)
    W <- make_unif_W(p, n.h, a)
    out_layer <- logistic(cbind(1, X) %*% W)
    return(out_layer)
}