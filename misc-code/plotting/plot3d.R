#' Plots only the cytogram (\code{ylist} with masses \code{countlist}) and not
##' the model.
##'
##' @param ylist List of covariates.
##' @param countslist List of counts. Defaults to null.
##' @param tt Time point.
##' @param ... Other arguments to plot3d.covarem().
##'
plot3d_cytogram <- function(ylist, X, countslist=NULL, tt, ...){
  plot3d.covarem(NULL, ylist, X, countslist, tt, show.xb.constraint = FALSE, ...)
}

##' Main 3d plotting function. It's possible to call this function without an
##' object |obj|. IN that case, it is a pure data plotting function.
##'
##' *Don't* do this.  !
##'
##' @param obj |covarem| object.
##' @param ylist List of covariates.
##' @param countslist List of counts. Defaults to null.
##' @param tt time point of interest.
##' @param show.xb.constraint If TRUE, show the ball constraint boundaries.
##'
##' @export
plot3d.covarem <- function(obj,
                           ## Understandably, data (ylist) might not be in the object.
                           ylist, countslist = NULL,
                           ## The time point of interest, out of 1:TT
                           tt,
                           ## Other options.
                           ## 2d scatterplot options
                           show.xb.constraint = FALSE,
                           cex.fac.2d = 1,
                           par_cex_2d = 1,
                           pt_col = rgb(0 ,0, 1, 0.1),
                           ## 3d scatterplot options
                           cex.fac.3d = 1,
                           ## 3d scatterplot options
                           destin = NULL
                           ){

  ## Define layout
  m = matrix(c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5,
               1, 1, 2, 2, 3, 3, 4, 4, 5, 5,
               6, 6, 6, 7, 7, 7, 7, 9, 9, 9,
               6, 6, 6, 8, 8, 8, 8, 9, 9, 9),
               nrow = 4, ncol = 10, byrow=TRUE)
  if(is.null(obj)){ m = m[1:2,] } ## Handling for missing obj; only the
                                  ## cytograms are to be plotted.
  layout(m)
  par(oma=c(3,1,2,1)) ## setting outer margin

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
    one_dim_scatterplot(ylist, obj, tt,
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
  ## Plot pies over time ###
  ##########################
  ## Todo: replace the following code with with plot_pie()
  pies = lapply(1:numclust, function(iclust){ obj$pie[,iclust] })
  pies.right.now = sapply(1:numclust, function(iclust){pies[[iclust]][[tt]]})
  names(pies.right.now) = paste("Clust", 1:numclust)
  cols = RColorBrewer::brewer.pal(numclust, "Set3")
  plot(NA, xlim=c(0,TT), ylim=c(0,1),
       ylab = "",
       xlab="",##"time, t=1,..,T",
       cex.axis = 2,
       axes=FALSE)
  add_date_ticks(obj)
  for(iclust in 1:numclust){
    lines(pies[[iclust]], col=cols[iclust], lwd=2)##, lty=iclust)
  }
  for(iclust in 1:numclust){
    x = TT-2*numclust+2*iclust
    text(x=x,
         y=pies[[iclust]][x],
         label = iclust, cex=1.5)
  }
  abline(v = tt, col='green', lwd=3)
  ## text(x=tt, y=1, label="Now", cex=3)
  abline(h=seq(from=0, to=1, by= 0.1), col='grey80')

  legend("topleft", legend="Cluster Probabilities", cex=3,
         bty = "n")

  ## TODO: try this instead.
  ## plot_pie(res)


  ####################
  ## Plot pies now ###
  ####################
  ## cols = rep("lightblue", numclust)
  names(pies.right.now) = 1:numclust
  barplot(pies.right.now, ylim=c(0,1),
          cex.axis = 2,
          cex.names = 2.5)
  abline(h=seq(from=0, to=1, by= 0.1), col='grey80')
  barplot(pies.right.now,
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
  ##   obj$loglikelihoods = objective(obj$mn, obj$pie, obj$sigma, ylist,
  ##                                  pie_lambda = obj$pie_lambda,
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
  ## mtext(text=paste0("Outer line "),side=1,outer=TRUE)

}


##' Makes a GENERIC 'fancy' plotly plot from algorithm output, for any two dimensions.
##' @param res result of running covarem().
##' @return plotly object.
fancyplot <- function(res, saveplot = FALSE, filename = NULL,
                      title = "3D Scatter plot",
                      dims = 1:2
                      ){

  ## Setup
  dimdat = res$dimdat
  assertthat::assert_that(all(dims %in% 1:dimdat))
  assertthat::assert_that(res$numclust == 4)

  ## Create a plottable data object
  data = lapply(1:res$numclust, function(iclust){
    x = res$mn[,dims[1],iclust]
    y = res$mn[,dims[2],iclust]
    z = 1:(res$TT)
    c = rep(iclust, res$TT)
    s = res$pie[,iclust]*50
    dat = data.frame(x, y, z, c, s)
    names(dat) = paste0( c("x", "y", "z", "c", "s"), iclust)
    dat
  })
  data = do.call(cbind, data)

  ## Setup
  scene = list(camera = list(eye = list(x = -1.25, y = 1.25, z = 1.25)),
               zaxis = list(title = "Time"),
               xaxis = list(title = "Dim 1"),
               yaxis = list(title = "Dim 2"))

  ## Make plot device
  p <- plotly::plot_ly(data, x = ~x1, y = ~y1, z = ~z1, type = 'scatter3d',
                       mode = 'lines+markers',
                       line = list(width = 2, fill = ~c1, colorscale = 'Viridis'),
                       marker = list(size = ~s1, fill = ~c1,
                                     colorscale = 'Viridis'),
                       name="Cluster 1")  %>%
  plotly::add_trace(x = ~x2, y = ~y2, z = ~z2,
            line = list(width = 2, fill = ~c2, colorscale = 'Viridis'),
            marker = list(size=~s2, fill=~c2),
            name = "Cluster 2")%>% ## ,
  plotly::add_trace(x = ~x3, y = ~y3, z = ~z3,
                       line = list(width = 2, fill = ~c3, colorscale = 'Viridis'),
                       marker = list(size = ~s3, fill = ~c3),
                       name = "Cluster 3")  %>%
  plotly::add_trace(x = ~x4, y = ~y4, z = ~z4,
                       line = list(width = 2, fill = ~c4, colorscale = 'Viridis'),
                       marker  =  list(size = ~s4, fill = ~c4),
                       name = "Cluster 4") %>%
  ## {if(has_name(data, "x5")) plotly::add_trace(x = ~x5, y = ~y5, z = ~z5,
  ##                      line = list(width = 2, fill = ~c5, colorscale = 'Viridis'),
  ##                      marker = list(size=~s5, fill=~c5),
  ##                      name="Cluster 5") else .} %>%
  ## {if(has_name(data, "x6")) plotly::add_trace(x = ~x6, y = ~y6, z = ~z6,
  ##                      line = list(width = 2, fill = ~c6, colorscale = 'Viridis'),
  ##                      marker = list(size=~s6, fill=~c6),
  ##                      name="Cluster 6") else .} %>%
  ## {if(has_name(data, "x7")) plotly::add_trace(x = ~x7, y = ~y7, z = ~z7,
  ##                      line = list(width = 2, fill = ~c7, colorscale = 'Viridis'),
  ##                      marker = list(size=~s7, fill=~c7),
  ##                      name="Cluster 7") else .} %>%
  ## {if(has_name(data, "x8")) plotly::add_trace(x = ~x8, y = ~y8, z = ~z8,
  ##                      line = list(width = 2, fill = ~c8, colorscale = 'Viridis'),
  ##                      marker = list(size=~s8, fill=~c8),
  ##                      name="Cluster 8") else .} %>%
  plotly::layout(title = title,
         scene = scene)
         ## autosize = F, width = 500, height = 500)
  if(saveplot) htmlwidgets::saveWidget(as.widget(p), filename)
  return(p)
}

##' Makes a table of coefficient values from an covarem class object |res|.
fancytable <- function(res, type=c("alpha", "beta")){
  if(type=="alpha")return(plotmat(tables$betas))
  if(type=="beta")return(plotmat(tables$alphas))

}

##' Helper to plot a table.
plotmat <- function(mat){
  p <- plot_ly(
    type = 'table',
    header = list(
        values = c("Dim", colnames(mat)),
        line = list(width = 1, color = 'black'),
        fill = list(color = 'rgb(235, 100, 230)'),
        font = list(family = "Arial", size = 14, color = "white")
    ),
    cells = list(
      values = rbind(
        rownames(mat),
        t(as.matrix(unname(mat)))
      ),
      align = c('left', rep('center', ncol(mat))),
      line = list(color = "black", width = 1),
      fill = list(color = c('rgb(235, 193, 238)', 'rgba(228, 222, 249, 0.65)')),
      font = list(family = "Arial", size = 12, color = c("black"))
    ))
  return(p)
}

##' Specific helper to obtain the coefficient table.
get_table <- function(res){
  ## Add table of betas
  betas = round(do.call(cbind, res$beta),2)
  rownames(betas) = paste("beta:", c("intrcpt", paste0("coef", 1:5)))
  betas[,seq(from = 1, to = ncol(betas), by=2)]
  betas = do.call(cbind, lapply(res$beta, function(mybeta){
    vec = rep(NA, 2*nrow(mybeta))
    vec[seq(from = 1, to = nrow(mybeta)*2, by=2)] = mybeta[,1]
    vec[seq(from = 2, to = nrow(mybeta)*2, by=2)] = mybeta[,2]
    return(vec)
  }))
  rownames(betas) = c("Intercept", "",
                     "SST", "",
                      "Salinity", "",
                     "Iron", "",
                     "Phosphorus", "",
                     "Chlorophyll", "")
  betas = signif(betas, 2)
  betas = cbind(rep(c(1,2),6), betas)
  colnames(betas) = c("Dim", paste0("Clust ", c(1,2,3,4)))

  ## Add table of alphas
  alphas = as.matrix(t(res$alpha))
  alphas = round(alphas,2)
  abg <- matrix(c(rep("grey90",4),
                  rep("grey80",4)), nrow=6, ncol=4,
                byrow=TRUE)
  rownames(alphas) = c("Intercept",
                     "SST",
                      "Salinity",
                     "Iron",
                     "Phosphorus",
                     "Chlorophyll")
  colnames(alphas) =  paste0("Clust ", c(1,2,3,4))
  return(list(betas=betas, alphas=alphas))
}

##' Use plotly to draw 3d figures of cytograms.
##'
##' (TODO: https://laustep.github.io/stlahblog/posts/plotly_trisurf.html : This
##' is a reference to try to make the ellipsoids actually elliptical)
##'
##' @param res covarem object.
##'
##' @return A plotly object.
plot3d_plotly.covarem <- function(res = NULL,
                               ylist,
                               countslist = NULL,
                               tt = NULL,
                               date = NULL,
                               fac = 2,
                               eye = list(x = 1.25, y = 1.25, z = 1.25),
                               xmax = NULL, ymax = NULL, zmax = NULL){

  numclust = res$numclust

  ################################################
  ## Make a matrix |dat| to use for plotting #####
  ################################################
  if(is.null(tt)){
    assertthat::assert_that(length(tt))
    tt = which(names(ylist) == date)
    assertthat::assert_that(length(tt) == 1)
  }
  y = ylist[[tt]]
  counts = countslist[[tt]]
  if(is.null(countslist)) counts = rep(1, nrow(y))
  stopifnot(nrow(y) == length(counts))
  dat = data.frame(cbind(y, counts))
  colnames(dat) = c("fsc_small", "chl_small", "pe", "count")
  if(!is.null(countslist)){
    dat[,"logcount"] = log(dat[,"count"] + 1) * fac
  } else {
    dat[,"logcount"] = 3
  }

  ## Set the axis limits
  if(is.null(xmax)) xmax = max(dat[,"fsc_small"])
  if(is.null(ymax)) ymax = max(dat[,"chl_small"])
  if(is.null(zmax)) zmax = max(dat[,"pe"])

  ############################
  ## Make base plot object ###
  ############################
  p = plotly::plot_ly() %>%
    ## Adding data points.
    plotly::add_trace(data = dat, type = 'scatter3d',
              x = ~fsc_small,
              y = ~chl_small,
              z = ~pe,
              mode = 'markers',
              marker = list(size = ~logcount,
                            color = 'rgba(0, 0, 253, 0.5)',
                            line = list(color = 'rgba(0, 0, 253, 0.5)',
                                        width = 0))) %>%
    ## Default Layout
    plotly::layout(scene = list(aspectmode = "data",
                                xaxis = list(range = c(0,xmax)),
                                yaxis = list(range = c(0,ymax)),
                                zaxis = list(range = c(0,zmax)),
                                camera = list(eye = eye)
                                ))

  ## ################################
  ## ## Add covariance information ##
  ## ################################
  ## if(!is.null(res)){
  ##   for(iclust in 1:numclust){
  ##       ellipse = rgl::ellipse3d(res$sigma[iclust,,],
  ##                                   centre = res$mn[tt,,iclust])
  ##       p = p %>% plotly::add_trace(size = 1,
  ##                                   x = ellipse$vb[1,], y = ellipse$vb[2,], z = ellipse$vb[3,],
  ##                                   opacity = 0.2, alphahull = 0,
  ##                                   type='mesh3d')
  ##   }
  ## }
  return(p)
}





##' Function generic for the main 3d plotting functionality.
plot3d <- function (x, ...) {
   UseMethod("plot3d", x)
}

##' Function generic for the plotly-type 3d plotting functionality.
plot3d_new <- function (x, ...) {
   UseMethod("plot3d_new", x)
}

##'
##' Helper to delete zero rows in a matrix.
##' @param mat Numeric matrix.
##'
##' @return Matrix with all zero rows removed.
reformat <- function(mat){
  inds = which(apply(mat, 1, function(myrow){
    all(myrow==0)
  }))
  if(any(inds)) return(mat[-inds,,drop=FALSE]) else return(mat)
}


##' Helper: to get ranges for each dimension, of the list of cytograms |ylist|.
##'
##' @param ylist List of cytogram data.
##'
##' @return Matrix whose rows contain the range (min & max) of each dimension.
get_range_from_ylist <- function(ylist){
  maxs = do.call(rbind, lapply(ylist, function(a) apply(a,2,max)))
  mins = do.call(rbind, lapply(ylist, function(a) apply(a,2,min)))
  maxs = apply(maxs, 2, max)
  mins = apply(mins, 2, min)
  return(rbind(mins, maxs))
}


##' Making *heatmap* for two dimensions of the original three (those in |ind|).
one_dim_heatmap <- function(ylist, obj, tt, countslist = NULL, dims = c(1,2)){

  ## Extract data
  y = ylist[[tt]][,dims]
  labs = colnames(ylist[[1]])
  maxcount = max(unlist(countslist))
  if(!is.null(obj)){
    mns = obj$mn
    TT = obj$TT
    numclust = obj$numclust
  }

  ## Get plot ranges
  ranges = get_range_from_ylist(ylist)
  ylim = ranges[,dims[2]]
  xlim = ranges[,dims[1]]##range(all.y[,dims[1]])
  ylab = labs[dims[2]]
  xlab = labs[dims[1]]

  ## Making the background color different
  allcounts = (sapply(countslist, sum))
  prop = allcounts[[tt]] / max(allcounts)
  yrange = ylim[2] - ylim[1]
  ylim_other = ylim[1]
  ylim[1] = ylim[1] - 1/5 * yrange

  ## Create empty plot
  main0 = ""
  plot(NA, ylim = ylim, xlim = xlim, cex = 3,
       ylab = ylab, xlab = xlab,
       cex.lab = 2,
       cex.axis = 2)
  title(main = main0, cex.main = 3)

  ## Add heatplot, using the *sparse* (full d^3) counts
  ymat = y_to_ymat(y, grid, dims)
  drawmat_precise2(ymat)

  make_ymat <- function(y, counts, grid, dims){

    ## Extract only the two columns of interest
    y = y[,dims]

    ## Then, collapse the data
    y = ybin_list[[100]]
    z = y[,1:2]
    dim(z)
    dim(unique(z))

    ## Make into a sparse grid.




    ## Then, collapse

  }

  ##  Input: a d x d matrix

  ## 1. Collapse the ylist.

  ylist[]
}

##' Making data plot for two dimensions of the original three (those in |ind|).
one_dim_scatterplot <- function(ylist, obj, tt, countslist = NULL, dims = c(1,2),
                                cex_fac = 10,
                                xlab = NULL,
                                ylab = NULL,
                                xlim = NULL,
                                ylim = NULL,
                                constant_total_count = FALSE,
                                ## The remaining arguments are for \code{one_dim_scatterplot_addmodel()}
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
    one_dim_scatterplot_addmodel(obj, tt, dims,
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

##' Helper function for \code{one_dim_scatterplot()}, to add means and ellipses
##' for confidence regions.
##'
##' @param obj A |covarem| class object.
##' @param tt time point.
##' @param dims e.g. \code{c(1,2)}, containing the two dimensions you want.
##' @param col color of the added means and confidence regions.
##'
##' @return NULL.
one_dim_scatterplot_addmodel <- function(obj, tt, dims,
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
  pies = lapply(1:numclust, function(iclust){ obj$pie[,iclust] })
  pies.right.now = sapply(1:numclust, function(iclust){pies[[iclust]][[tt]]})
  mn.cex = pies.right.now/max(pies.right.now) * mn_cex_fac

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

  ## Collect pies
  if(!is.null(obj)){
    numclust = obj$numclust
    pies = lapply(1:numclust, function(iclust){ obj$pie[,iclust] })
    pies.right.now = sapply(1:numclust, function(iclust){pies[[iclust]][[tt]]})
    mn.cex = pies.right.now/max(pies.right.now)*8
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
