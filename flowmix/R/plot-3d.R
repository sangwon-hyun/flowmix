##' Main 3d plotting function. It's possible to call this function without an
##' object |obj|. IN that case, it is a pure data plotting function.
##'
##' *Don't* do this.
##'
##' @param obj flowmix object.
##' @param ylist List of covariates.
##' @param countslist List of counts. Defaults to null.
##' @param tt time point of interest.
##' @param show.xb.constraint If TRUE, show the ball constraint boundaries.
##'
##' @export
plot_3d <- function(obj,
                    ylist, countslist = NULL, ## The time point of interest, out of 1:TT
                    tt, ## Other options.
                    ## 2d scatterplot options
                    show.xb.constraint = FALSE, cex.fac.2d = 1, par_cex_2d = 1,
                    pt_col = rgb(0 ,0, 1, 0.1), ## 3d scatterplot options
                    cex.fac.3d = 1, ## 3d scatterplot options
                    destin = NULL){

  ## Define layout
  m = matrix(c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5,
               1, 1, 2, 2, 3, 3, 4, 4, 5, 5,
               6, 6, 6, 7, 7, 7, 7, 9, 9, 9,
               6, 6, 6, 8, 8, 8, 8, 9, 9, 9),
               nrow = 4, ncol = 10, byrow=TRUE)
  if(is.null(obj)){ m = m[1:2,] } ## Handling for missing obj; only the
                                  ## cytograms are to be plotted.
  layout(m)
  par(oma = c(3, 1, 2, 1)) ## setting outer margin

  ## Setup
  TT = length(ylist)
  assert_that(tt %in% 1:TT)
  all.y = do.call(rbind, ylist)
  only_plot_cytograms = (is.null(obj))
  if(!only_plot_cytograms){
    obj = reorder_clust(obj)
    mns = obj$mn
    numclust = obj$numclust
    p = ncol(obj$X)
  }

  ## Scale the biomass (|countslist|) by the total biomass in that cytogram.
  counts_sum = sapply(countslist, sum)
  fac = median(counts_sum)
  countslist = lapply(countslist, function(counts)counts/sum(counts) * fac)


  ###############################
  ## Make the three data plots ##
  ###############################
  ## par(mar=c(1,1,3,1))
  par(mar = c(5.1, 5.1, 4.1, 2.1))
  dimslist = list(1:2, 2:3, c(3,1))
  for(dims in dimslist){
    ## one_dim_scatterplot(ylist, obj, tt,
    scatterplot_2d(ylist, obj, tt,
                   countslist = countslist,
                   dims = dims,
                   cex_fac = cex.fac.2d,
                   pt_col = pt_col,
                   lwd = 2)
  }
  par(cex=0.8)

  par(mar=c(1,1,3,1))
  phis = c(10,50)
  for(phi in phis){
    one_3d_plot(ylist, obj, tt, countslist = countslist, phi = phi,
                cex.fac = cex.fac.3d)
  }


  ## Add map with cruise location.
  make_map(obj, tt, destin = destin)
  ## make_map_ggplot(obj, tt, destin = destin)

  ######################
  ## Plot covariates ###
  ######################
  ## plot(NA, xlim=c(0,TT*1), ylim = ylim.cov, ylab = "Covariates", xlab="",
  ##      axes = FALSE)
  plot_covariates(obj, tt=tt)

  ## If there is no |obj|, stop here.
  if(only_plot_cytograms){ return(NULL) }

  ##########################
  ## Plot probs over time ###
  ##########################
  probs = lapply(1:numclust, function(iclust){ obj$prob[,iclust] })
  probs.right.now = sapply(1:numclust, function(iclust){probs[[iclust]][[tt]]})
  names(probs.right.now) = paste("Clust", 1:numclust)
  cols = RColorBrewer::brewer.pal(numclust, "Set3")
  plot(NA, xlim=c(0,TT), ylim=c(0,1),
       ylab = "",
       xlab="",##"time, t=1,..,T",
       cex.axis = 2,
       axes=FALSE)
  add_date_ticks(obj)
  for(iclust in 1:numclust){
    lines(probs[[iclust]], col=cols[iclust], lwd=2)##, lty=iclust)
  }
  for(iclust in 1:numclust){
    x = TT-2*numclust+2*iclust
    text(x=x,
         y=probs[[iclust]][x],
         label = iclust, cex=1.5)
  }
  abline(v = tt, col='green', lwd=3)
  ## text(x=tt, y=1, label="Now", cex=3)
  abline(h=seq(from=0, to=1, by= 0.1), col='grey80')

  legend("topleft", legend="Cluster Probabilities", cex=3,
         bty = "n")

  ####################
  ## Plot probs now ###
  ####################
  ## cols = rep("lightblue", numclust)
  names(probs.right.now) = 1:numclust
  barplot(probs.right.now, ylim=c(0,1),
          cex.axis = 2,
          cex.names = 2.5)
  abline(h=seq(from=0, to=1, by= 0.1), col='grey80')
  barplot(probs.right.now,
          col = cols,
          lwd = 1, add=TRUE,
          names.arg = rep("", numclust))
  legend("topright", legend=paste0("Cluster Probabilities"), cex=3,
         bty = "n")
  legend("topright", legend=paste("\nTime",tt, "out of", TT), cex=3,
         bty = "n")



  ## ## Add total biomass (sum)
  ## plot(counts_sum/max(counts_sum), lwd=3, axes=FALSE, ylab="", xlab="", type='l')
  ## title(main="Cytogram total biomass over time", cex.main=2)
  ## ## text(y=(counts_sum/max(counts_sum))[1], x=10,
  ## ##      label="total\nbiomass")
  ## add_date_ticks(obj)
  ## abline(v = tt, col = 'green')

  ## ###########################################
  ## ## Add the likelihood of each time point ##
  ## ###########################################
  ## if(is.null(obj$loglikelihoods)){
  ##   ## Also get per-cytogram likelihoods
  ##   obj$loglikelihoods = objective(obj$mn, obj$prob, obj$sigma, ylist,
  ##                                  prob_lambda = obj$prob_lambda,
  ##                                  mean_lambda = obj$mean_lambda,
  ##                                  alpha = obj$alpha, beta = obj$beta,
  ##                                  countslist = obj$countslist,
  ##                                  each = TRUE)
  ## }
  ##   ## ntlist = sapply(ylist, length)
  ##   if(is.null(countslist)){
  ##     ntlist = sapply(countslist, sum)
  ##   } else {
  ##     ntlist = sapply(ylist, nrow)
  ##   }
  ## plot(obj$loglikelihoods / ntlist, type='o', axes=FALSE,
  ##      cex.lab = 2, ylab = "In-sample Log likelihoods",
  ##      xlab = "time, t=1,..,T")
  ## axis(1); axis(2)
  ## title(main="Likelihoods for each cytogram", cex.main = 2)
  ## abline(v = tt, col = 'green')


  ## Add outer title
  main0 = paste0("Time ", tt, " out of ", TT)
  mtext(bquote(bold(.(main0))), line=-1, side=3, outer=TRUE, cex=3,
        font.main = 1)
  time_now = rownames(obj$X)
  time_now = gsub("T", " ", time_now)
  mtext(bquote(bold(.(time_now))), line=-4, side=3, outer=TRUE, cex=3,
        font.main = 1)
}



##' Function generic for the main 3d plotting functionality.
plot3d <- function (x, ...) {
   UseMethod("plot3d", x)
}


##' Making data plot for two dimensions of the original three (those in |ind|).
##'
##' @export
scatterplot_2d <- function(ylist, obj, tt, countslist = NULL, dims = c(1,2),
                           cex_fac = 10,
                           xlab = NULL,
                           ylab = NULL,
                           xlim = NULL,
                           ylim = NULL,
                           constant_total_count = FALSE,
                           ## The remaining arguments are for \code{scatterplot_2d_addmodel()}
                           add_clust_labels = TRUE,
                           only_model = FALSE,## If TRUE, only draw the model.
                           subset_of_clust = NULL, ## if not null, contains the clusters to plot
                           col = NULL,
                           lty = 2,
                           lwd = 1/2,
                           not_so_small_cex = FALSE,
                           mn_cex_fac = 5,
                           pt_col = rgb(0 ,0, 1, 0.1),
                           par_cex = 1
                           ){

  ## Extract data
  par(cex=par_cex)
  y = ylist[[tt]][,dims]
  labs = colnames(ylist[[1]])
  if(!is.null(countslist)){
    maxcount = max(unlist(countslist))
    ## countslist = lapply(countslist, function(counts){
    ##   counts / maxcount
    ## })
    counts = countslist[[tt]] /maxcount
    if(constant_total_count) counts = counts/sum(counts)
  } else {
    ## countslist = lapply(ylist, function(y){ rep(1, nrow(y)) })
    counts = rep(1, nrow(y))
  }

  ## Aggregating counts into the two dimensions
  yy = collapse_3d_to_2d(ylist[[tt]], counts, dims)
  counts = yy[,3]
  y = yy[,1:2]

  if(!is.null(obj)){
    mns = obj$mn
    TT = obj$TT
    numclust = obj$numclust
  }

  ## Get plot ranges
  ranges = get_range_from_ylist(ylist)
  if(is.null(ylim)) ylim = ranges[,dims[2]]
  if(is.null(xlim)) xlim = ranges[,dims[1]]##range(all.y[,dims[1]])
  if(is.null(ylab))ylab = labs[dims[2]]
  if(is.null(xlab))xlab = labs[dims[1]]

  ## Make the 2d plot
  ## if(!is.null(countslist)){
  ##   cex = sqrt(count) / max_sqrt_count * cex.fac
  ## }

  ## browser()
  ## ifelse(is.null(countslist), NULL, counts)
  ## if(is.null(countslist)) counts = NULL

  one_2d_plot(y = y,
              counts = counts,
              xlim = xlim, ylim = ylim,## cex=cex,
              ylab = ylab, xlab = xlab,
              pt_col = pt_col,
              cex_fac = cex_fac)

  ## Add the means
  if(!is.null(obj)){
    scatterplot_2d_addmodel(obj, tt, dims,
                                 add_clust_labels = add_clust_labels,
                                 subset_of_clust = subset_of_clust, col=col,
                                 lwd = lwd,
                                 lty = lty,
                                 not_so_small_cex = not_so_small_cex,
                                 mn_cex_fac = mn_cex_fac
                                 )
  }

  par(cex=1)
}

##' Helper function for \code{scatterplot_2d()}, to add means and ellipses
##' for confidence regions.
##'
##' @param obj A |covarem| class object.
##' @param tt time point.
##' @param dims e.g. \code{c(1,2)}, containing the two dimensions you want.
##' @param col color of the added means and confidence regions.
##'
##' @return NULL.
scatterplot_2d_addmodel <- function(obj, tt, dims,
                                         col = "tomato",
                                         add_clust_labels = TRUE,
                                         clust_labels = NULL,
                                         ellipse = TRUE,
                                         subset_of_clust = NULL,
                                         lwd = 1/2,
                                         lty = 2,
                                         not_so_small_cex = FALSE,
                                         mn_cex_fac = 5
                                         ){

  ## Take a few objects
  mns = obj$mn
  TT = obj$TT
  numclust = obj$numclust

  ## Basic checks
  if(is.null(col)) col = "tomato"
  stopifnot(all(dims%in%c(1:obj$dimdat)))
  stopifnot(length(dims)==2)
  if(is.null(subset_of_clust)){
    allclusts = clust_labels = 1:numclust
  } else {
    stopifnot(all(subset_of_clust %in% 1:numclust))
    allclusts = subset_of_clust
  }
  if(add_clust_labels){
    if(!is.null(clust_labels)){
      stopifnot(length(clust_labels)==length(allclusts))
    }
      ## clust_labels = 1:numclust
  }

  ## Add the ellipses for the covariances
  if(ellipse){
  for(iclust in allclusts){
    lines(ellipse::ellipse(x = obj$sigma[iclust,dims,dims],
                           centre = mns[tt,dims, iclust]),
          lwd = lwd,
          col = col,
          lty = lty)
  }
  }

  ## Add fitted means
  probs = lapply(1:numclust, function(iclust){ obj$prob[,iclust] })
  probs.right.now = sapply(1:numclust, function(iclust){probs[[iclust]][[tt]]})
  ### TODO: probs.right.now = sqrt(probs.right.now)
  mn.cex = probs.right.now/max(probs.right.now) * mn_cex_fac

  ## Make the minimum size not so small
  if(not_so_small_cex){
    mn.cex = pmax(mn.cex, 2 * mn_cex_fac/5)
  }

  for(ii in 1:length(allclusts)){
  ## for(iclust in allclusts){
    iclust = allclusts[ii]


    ## Add means
    points(x = mns[tt,dims[1],iclust],
           y = mns[tt,dims[2],iclust],
           col = col,
           pch=16, cex=mn.cex[iclust])
  }

  for(ii in 1:length(allclusts)){

    iclust = allclusts[ii]

    ## Add cluster labels on means
    if(add_clust_labels){
      if(!is.null(allclusts)){
        label = clust_labels[ii]
      } else {
        label = iclust
      }
      text(x = mns[tt,dims[1],iclust],
           y = mns[tt,dims[2],iclust],
           labels = label,
           col = 'black', pch = 16, cex = 2.5 * sqrt(mn_cex_fac) / sqrt(5),
           font = 2, pos = 4)
    }
  }

}



##' Making a 3d scatter plot with a certain angle..
##' @param ... Additional arguments to plot3d::scatter3d().
one_3d_plot <- function(ylist, obj=NULL, tt, countslist=NULL, phi = 40,
                        cex.fac = 1,
                        cex.axis = 1,
                        xlab = NULL,
                        ylab = NULL,
                        zlab = NULL,
                        ticktype = "detailed",
                        col = rgb(0,0,1,0.1),
                        ...){

  y = ylist[[tt]]
  if(is.null(xlab)) xlab = colnames(y)[1]
  if(is.null(ylab)) ylab = colnames(y)[2]
  if(is.null(zlab)) zlab = colnames(y)[3]
  if(is.null(countslist)){
    cex = 1 * cex.fac
  } else {
    maxcount = max(unlist(countslist))
    cex = countslist[[tt]]
    cex = cex / maxcount
    cex = cex * cex.fac
    ## cols = sapply(cex, function(cx) rgb(0, 0, 1, cx))
  }


  ranges = get_range_from_ylist(ylist)
  xlim = ranges[,1]##range(all.y[,dims[1]])
  ylim = ranges[,2]
  zlim = ranges[,3]

  ## Collect probs
  if(!is.null(obj)){
    numclust = obj$numclust
    probs = lapply(1:numclust, function(iclust){ obj$prob[,iclust] })
    probs.right.now = sapply(1:numclust, function(iclust){probs[[iclust]][[tt]]})
    mn.cex = probs.right.now/max(probs.right.now)*8
    mns = obj$mn
  }


  ## Using scatter3D package, make plot:


  ## Plot the data points.
  plot3D::scatter3D(y[,1], y[,2], y[,3],
            col = col,
            pch = 16,
            bty='g', phi=phi,
            xlim = xlim, ylim = ylim, zlim = zlim,
            xlab = xlab, ylab = ylab, zlab = zlab,
            cex = sqrt(cex) * 10,
            cex.lab = 2, ## Trying to get the labels to magnify
            ticktype = ticktype,## Using detailed ticks
            cex.axis = cex.axis,
            ...) 

  if(!is.null(obj)){

  ## Add labels.
  plot3D::text3D(mns[tt,1,],
                 mns[tt,2,],
                 mns[tt,3,],
                 labels = 1:numclust,
                 add = TRUE,
                 colkey = FALSE,
                 cex = 2.5)

  ## Add cluster centers.
  mn.pch = 16##"O"
  plot3D::scatter3D(mns[tt,1,],
                    mns[tt,2,],
                    mns[tt,3,],
                    col = "tomato",##"orange", ## "red",
                    pch = mn.pch,
                    add = TRUE,
                    cex = mn.cex,## * 0.5,
                    type = 'h',
                    lwd = 3)##, type="h", pch=16)
  }
}




##' A simple function for plotting 3-dimensional data.
##'
plot3d_simple <- function(obj,
                          ## Understandably, data (ylist) might not be in the object.
                          ylist, countslist = NULL,
                          ## The time point of interest, out of 1:TT
                          tt,
                          ## Other options.
                          ## 2d scatterplot options
                          cex.fac.2d = 1,
                          pt_col = rgb(0 ,0, 1, 0.1),
                          ## 3d scatterplot options
                          cex.fac.3d = 1
                          ){

  ## Define layout
  par(mfrow = c(1,3))
  par(oma = c(3, 1, 2, 1)) ## setting outer margin

  ## Setup
  TT = length(ylist)
  assert_that(tt %in% 1:TT)
  all.y = do.call(rbind, ylist)
  only_plot_cytograms = (is.null(obj))
  if(!only_plot_cytograms){
    obj = reorder_clust(obj)
    mns = obj$mn
    numclust = obj$numclust
    p = ncol(obj$X)
  }

  ## Temporary
  dimnames = colnames(ylist[[1]])

  ## Scale the biomass (|countslist|) by the total biomass in that cytogram.
  counts_sum = sapply(countslist, sum)
  fac = median(counts_sum)
  countslist = lapply(countslist, function(counts)counts/sum(counts) * fac)


  ###############################
  ## Make the three data plots ##
  ###############################
  ## par(mar=c(1,1,3,1))
  par(mar = c(5.1, 5.1, 4.1, 2.1))
  dimslist = list(1:2, 2:3, c(3,1))
  for(dims in dimslist){
    scatterplot_2d(ylist = ylist,
                   countslist = countslist,
                   obj = obj,
                   tt = tt,
                   dims = dims, cex_fac=cex.fac.2d,
                   pt_col = rgb(0 ,0, 1, 0.1),
                   xlab = dimnames[dims[1]],
                   ylab = dimnames[dims[2]])
  }
}
