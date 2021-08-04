##' For 2-dimensional data, make a single scatterplot for one time points' worth
##' of data.
##'
##' CAN'T handle the binned case (with \code{countslist}) yet.
##'
##' @param ylist Data.
##' @param tt Time point.
##' @param res flowmix object.
##' @param drawslist A list containing the randomly generated memberships using
##'   the posterior cluster probabilities, which can be generated using
##'   \code{draw_membership()}. Defaults to NULL.
##' @param resp Posterior probabilities of cluster memberships
##'   (responsibilities). Defaults to NULL.
##' @param iclust Cluster number(s). Defaults to NULL.
##' @param cols Colors for each cluster. Defaults to NULL.
##'
##' @importFrom magrittr %>%
##'
##' @export
##'
plot_2d <- function(ylist, tt, res = NULL, drawslist = NULL, resp = NULL, iclust = NULL, cols = NULL){

  ## Setup
  stopifnot(ylist %>% .[[1]] %>% ncol() == 2)
  rngs = do.call(rbind, ylist) %>% apply(., 2, range)
  numclust = res$numclust
  if(is.null(iclust)){ numclusts = 1:numclust } else {numclusts = iclust}
  if(is.null(cols)) cols = RColorBrewer::brewer.pal(numclust, "Set3")[1:numclust]
  if(!is.null(cols)) assertthat::assert_that(length(cols) == numclust)
  if(!is.null(res)) assertthat::assert_that(is(res, "flowmix"))


  cartoon_points <- function(y, x, cex, col, pch = 16, ...){
    points(y = y, x = x, cex = cex * 10 + 2, pch = pch, col = "black", ...)
    points(y = y, x = x, cex = cex * 10 + 1, pch = pch, col = col)
  }


  ## Draw membership
  if(is.null(drawslist)){
    if(is.null(resp)) resp <- Estep(res$mn, res$sigma, res$prob, ylist = ylist,
                                     numclust = res$numclust, first_iter = TRUE)
    drawslist = draw_membership(resp)
  }
  draws = drawslist[[tt]]

  ## Make plot
  plot(NA, ylim = rngs[,2], xlim = rngs[,1], xlab="", ylab="")

  ## Add the data points
  for(iclust in numclusts) points(ylist[[tt]][which(draws[,iclust]==1),1:2],
                                  col =cols[iclust] %>% adjustcolor(alpha.f = 0.5),
                                  pch = 16)
  ## Add the cluster centers
  for(iclust in numclusts){
    cartoon_points(y = res$mn[tt,2,iclust],
                   x = res$mn[tt,1,iclust],
                   cex = res$prob[tt,iclust],
                   col = cols[iclust] %>% adjustcolor(alpha.f = 1))
    lines(ellipse::ellipse(x = res$sigma[iclust,,], centre = res$mn[tt,, iclust]),
          lwd = 2, col = cols[iclust], lty = 1)
  }
}



##' 2d Plotting functionality using ggplot2 (only 2d data).
##'
##' @param datobj_2d A data matrix with X on grid; three columns are assumed;
##'   the first two are the coordinates of the 2d data grid (e.g. "diam" and
##'   "chl"); the thrid column is named "counts".
##'
##' @return A ggplot object.
##'
##' @export
##'
bin_plot_2d <- function(datobj_2d, mn=NULL, sigma=NULL, prob=NULL, labels = NULL, fac = 20, colours = NULL){

  ## Basic checks
  if(is.null(colours)) colours = c("white", "blue")

  ## Get variable names
  varnames = datobj_2d %>% colnames()
  stopifnot(length(varnames) == 3) ## two data columns, one called "counts"
  stopifnot(varnames[3] == "counts")
  varname1 = varnames[1]
  varname2 = varnames[2]

  ## Information about clusters
  dt = data.frame(x = mn[1,], y = mn[2,], prob = prob * fac)
  if(!is.null(mn)){
    numclust = dim(mn)[2]
    if(is.null(labels)) labels = 1:numclust
    dt = cbind(dt, label = labels)
  }

  p = datobj_2d %>%
    ggplot() +
    theme_minimal() +
    geom_raster(aes(x = !!sym(varname1), y=!!sym(varname2), fill = counts)) +
    scale_fill_gradientn(colours = colours, guide="colorbar")+
    xlim(c(0,8)) + ylim(c(0, 8)) +
    theme(legend.position = "none")

  ## Add model.
  if(!is.null(mn)){
    for(iclust in 1:numclust){
      el = ellipse::ellipse(x = sigma[iclust,,], centre = mn[,iclust]) %>% as_tibble()
      p = p + geom_path(aes(x = x, y = y), data = el, colour = "red", lty = 2,
                        lwd = pmin(prob[iclust] * 5, 0.5))
    }

    ## Add points
    p = p + geom_point(aes(x = x, y = y, size = prob),
                       data = dt, colour = 'red') +
                    ## data = dt, colour = 'black', fill = 'red', shape = 21, stroke = 2) +
    scale_size_identity()

    ## Add labels
    p = p + geom_text(aes(x = x, y = y, label = label), size = 5, hjust = 0, vjust = 0, data = dt,
                      fontface = "bold", col='black')
  }

  return(p)
}



##' Plot three 2d panels of data, optionally with a model (from \code{obj}).
##'
##' @param obj flowmix object. If NULL, only data is drawn.
##' @param ylist Data.
##' @param countslist Defaults to NULL.
##' @param obj A flowmix object.
##' @param tt time point of interest, out of 1 through \code{length(ylist)}.
##'
##' @return A grob object containing a 3-panel plot.
##'
##' @export
plot_2d_threepanels <- function(obj = NULL, ## Understandably, data (ylist) might not be in the object.
                                ylist,
                                countslist = NULL, ## The time point of interest, out of 1:TT
                                tt,
                                labels = NULL,
                                plist_return = FALSE,
                                colours = NULL){

  ## Basic checks
  if(!is.null(obj)) stopifnot("flowmix" %in% class(obj))

  ## Setup
  TT = length(ylist)
  assertthat::assert_that(tt %in% 1:TT)
  mn = sigma = prob = NULL

  ## Scale the biomass (|countslist|) by the total biomass in that cytogram.
  counts_sum = sapply(countslist, sum)
  countslist = lapply(countslist, function(counts)counts/sum(counts))

  ###############################
  ## Make the three data plots ##
  ###############################
  dimslist = list(c(1:2), c(2:3), c(3,1))
  plist = lapply(dimslist, function(dims){

    ## Collapse to 2d
    y = ylist[[tt]]
    counts = countslist[[tt]]
    datobj_2d = collapse_3d_to_2d(y = y, counts = counts, dims = dims) %>% as_tibble()

    ## Extract model data
    ## obj = bestres
    if(!is.null(obj)){
      mn = obj$mn[tt,dims,]
      sigma = obj$sigma[,dims,dims]
      prob = obj$prob[tt,]
    }

    ## Make the heatmap
    p = bin_plot_2d(datobj_2d, mn, sigma, prob, labels, colours = colours)

    return(p)
  })

  if(plist_return) return(plist)

  ###############################
  ## Return a 3 x 1 data panel ##
  ###############################
  mydatetime = names(ylist)[tt]
  main_text_grob = grid::textGrob(mydatetime,
                                  gp = grid::gpar(fontsize = 20, font = 3))
  gridExtra::arrangeGrob(plist[[1]], plist[[2]], plist[[3]],
                          ncol = 3, top = main_text_grob)
}
