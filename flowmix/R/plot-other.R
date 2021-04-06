##' Plots the time series of covariates from an |obj| object.
##'
##' @param obj An object of class |flowmix|.
##' @param tt Time point to highlight.
##'
##' @return no return
##'
##' @export
plot_covariates <- function(obj, tt = NULL){

  ## Setup
  X = obj$X
  ylim.cov = range(X) * c(1, 0.7)
  cols = 1:5

  ## Produce plot
  par(mar = c(5,5,2,2),
      cex.axis = 2, cex.lab = 2)
  TT = nrow(X)
  matplot(X,
          col = 'grey',
          lwd = 0.5,
          type = 'l',
          xlim = c(0,TT* 1),
          ylim = ylim.cov,
          axes = FALSE,
          ylab = "",
          xlab = "")
  legend("topleft",
         legend = "Environmental Covariates",
         cex = 3,
         bty = "n")

  ## Make a green vertical line to signify the current time point.
  if(!is.null(tt)) abline(v = tt, col = 'green', lwd = 3)

  ## Add date ticks.
  if(!is.null(rownames(obj$X)) & lubridate::is.Date(rownames(obj$X)[1])){
    add_date_ticks(obj)
  } else {
    axis(2); axis(1);
  }

}


##' Plot cluster probabilities over time.
##'
##' @param res flowmix object.
##' @param iclusts Optionally, provide the cluster numbers to plot only a
##'   subset.
##'
##' @return NULL
##'
##' @export
plot_prob <- function(res, iclusts=NULL, main=NULL,
                     cols = NULL
                     ){

  ## Setup
  if(is.null(iclusts)) iclusts = c(1:res$numclust)
  if(is.null(cols)){ cols = RColorBrewer::brewer.pal(res$numclust, "Set3") }

  ## Reorder clusters in decreasing order of diam
  res <- reorder_clust(res)

  ## Make plot
  matplot(NA,
          xlim = c(0, res$TT),
          ylab = "",
          xlab = "",
          ylim = c(0, 1),
          axes = FALSE)
  abline(h = seq(from = 0, to = 1, by = 0.1), col='grey90', lwd=2, lty=3)
  matlines(res$prob[,iclusts], type = 'l', lty = 1, lwd = 3, col = cols[iclusts])

  ## Add a main title
  if(is.null(main)){
    title(main = "Cluster probabilities", cex.main=2)
  } else {
    title(main = main, cex.main=1)
  }

  ## Add date ticks
  if(!is.null(rownames(res$X)) & lubridate::is.Date(rownames(res$X)[1])){
    add_date_ticks(res)
  } else {
    axis(2);axis(1);
  }
}


##' For a matrix of CV scores (which are included in output from the function
##' \code{aggregateres()} or \code{blockcv_summary()}, make a 2d heatmap
##'
##' @param cvscore.mat Matrix containing CV scores.
##'
##' @export
plot_cvscore <- function(cvscore.mat){
  mat = cvscore.mat
  colnames(mat) = signif(as.numeric(colnames(mat)), 2)
  rownames(mat) = signif(as.numeric(rownames(mat)), 2)
  drawmat_precise_helper(mat, contour = FALSE,
                         ylab = expression(lambda[alpha]),
                         xlab = expression(lambda[beta]))
}


##' Another helper function to /precisely/ draw the entries of a matrix.
##' @param mat Matrix of interest.
##' @param contour If \code{TRUE}, draw a contour using
##'   \code{lattice::levelplot()}.
##' @param ... Other arguments to \code{lattice::levelplot()}.
##'
##' @return lattice object.
drawmat_precise_helper <- function(mat, contour = FALSE, ...){

## Dummy data
## data <- matrix(runif(100, 0, 5) , 10 , 10)

  if(is.null(colnames(mat))){
    colnames(mat) <- paste(rep("col\n",ncol(mat)),
                            c(1:ncol(mat)) , sep=" ")
    rownames(mat) <- paste(rep("row",nrow(mat)),
                           c(1:nrow(mat)) , sep=" ")
  }

  ## Color function
  colfun = colorRampPalette(c("blue", "red"))

  # plot it flipping the axis
  lattice::levelplot(t(mat[c(nrow(mat):1) , ]),
                     col.regions = colfun(100),
                     contour = contour,
                     ## xaxt = 'n',
                     las = 2,
                     ...)
}
