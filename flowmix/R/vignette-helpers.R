##' Add to 2d plot the estimated cluster mean and 95% confidence region, for
##' each cluster.
##'
##' @param p ggplot of 2d scatterplot of data
##' @param mn mean
##' @param sigma covariance
##' @param prob relative abundance
##' @export
##' @keywords internal
add_model_2d <- function(p, mn, sigma, prob){

  ## Basic checks
  stopifnot(nrow(mn) == 2)

  ## Harmlessly silencing some check() warnings
  x = NULL
  y = NULL
  label = NULL

  fac = 20
  dt = data.frame(x = mn[1,], y = mn[2,], prob = prob * fac)
  numclust = dim(mn)[2]
  labels = 1:numclust
  dt = cbind(dt, label = labels)


  for(iclust in 1:numclust){
    el = ellipse::ellipse(x = sigma[iclust,,], centre = mn[,iclust]) %>% as_tibble()
    p = p + geom_path(aes(x = x, y = y), data = el,
                      col = 'red',
                      lty = 2,
                      lwd = pmin(prob[iclust] * 8, 0.8))
  }

  ## Add points
  ## p = p + geom_point(aes(x = x, y = y, size = (prob)),
  ##                    data = dt, colour = 'red') +

  p = p + geom_point(aes(x = x, y = y, size = (prob)),
                     data = dt,
                     col = 'red') +
                     ## colour = mn_colours) +
    scale_size_area()


  ## Add labels
  ## p = p + geom_text(aes(x = x, y = y, label = label), size = 5, hjust = 0, vjust = 0, data = dt,
  ##                   fontface = "bold", col='black')

  ## Improve label style
  cex = rel(3)##ifelse(mn_colours==rgb(0,0,0,0.5), rel(3), rel(4))
  p = p + ggrepel::geom_text_repel(aes(x = x, y = y, label = label, point.size = sqrt(prob)),
                                   col = "red",
                                   ## col = mn_colours,
                                   cex = cex,
                                   bg.color = "white",
                                   bg.r = 0.1,
                                   fontface = "bold",
                                   ## point.size = NA,
                                   force_pull   = 5, # do not pull toward data points
                                   data = dt,
                                   seed = 1)
  p = p + theme(legend.position="none")
  return(p)
}




##' For a matrix of CV scores (which are included in output from the function
##' \code{aggregateres()} or \code{blockcv_summary()}, make a 2d heatmap
##'
##' @param cvscore.mat Matrix containing CV scores.
##'
##' @export
##' @keywords internal
plot_cvscore <- function(cvscore.mat){
  mat = cvscore.mat
  colnames(mat) = signif(as.numeric(colnames(mat)), 2)
  rownames(mat) = signif(as.numeric(rownames(mat)), 2)
  drawmat_precise_helper(mat, contour = FALSE,
                         ylab = expression(lambda[alpha]),
                         xlab = expression(lambda[beta]))
}


##' Another helper function to /precisely/ draw the entries of a matrix.
##'
##' @param mat Matrix of interest.
##' @param contour If \code{TRUE}, draw a contour using
##'   \code{lattice::levelplot()}.
##' @param ... Other arguments to \code{lattice::levelplot()}.
##'
##' @return lattice object.
##' @export
##' @keywords internal
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
  colfun = grDevices::colorRampPalette(c("blue", "red"))

  # plot it flipping the axis
  lattice::levelplot(t(mat[c(nrow(mat):1) , ]),
                     col.regions = colfun(100),
                     contour = contour,
                     ## xaxt = 'n',
                     las = 2,
                     ...)
}
