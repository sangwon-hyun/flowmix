##' Temporary function for plotting three different 3d models together. Useful
##' but not included in the package.
plot3d_compare.covarem <- function(obj1, obj2, obj3,
                                   ## Understandably, data (ylist) might not be in the object.
                                   ylist, countslist = NULL,
                                   ## The time point of interest, out of 1:TT
                                   tt,
                                   ## Other options.
                                   cex.fac = 1,
                                   destin = NULL ## Output of obj1. Used to
                                                   ## get the latlontime.Rdata
                                                   ## file, for plotting the map.
                                   ){

  ## Define layout
  m = matrix(c(1,2,3,4,5,
               6,7,8,9,10,
               11,12,13,14,15,
               16,17,17,17,17),
               nrow = 4, ncol = 5, byrow=TRUE)
  layout(m)
  par(oma=c(3,1,2,1)) ## setting outer margin

  ## Setup
  TT = length(ylist)
  assert_that(tt %in% 1:TT)
  all.y = do.call(rbind, ylist)
  obj1 = reorder_clust(obj1)
  obj2 = reorder_clust(obj2)
  obj3 = reorder_clust(obj3)

  ## Scale the biomass (|countslist|) by the total biomass in that cytogram.
  counts_sum = sapply(countslist, sum)
  fac = median(counts_sum)
  countslist = lapply(countslist, function(counts)counts/sum(counts) * fac)

  ###############################
  ## Make the 3 + 2 data plots ##
  ###############################
  for(jj in 1:3){
    obj = list(obj1, obj2, obj3)[[jj]]
    par(mar = c(5.1, 5.1, 4.1, 2.1))
    dimslist = list(1:2, 2:3, c(3,1))
    for(ii in 1:length(dimslist)){
      dims = dimslist[[ii]]
      one_dim_scatterplot(ylist, obj, tt,
                          countslist = countslist,
                          dims = dims,
                          cex.fac = cex.fac)
      if(ii==1){
        if(jj==1)legend("topleft", legend=paste0("Time ", tt, " out of ", TT), cex=2.5, bty = "n")
        legend("topleft", legend=paste0("\n", obj$numclust, " Clusters"), cex=2.5, bty = "n")
      }
    }
    par(mar=c(1,1,3,1))
    phis = c(10,50)
    for(phi in phis){
      one_3d_plot(ylist, obj, tt, countslist = countslist, phi = phi,
                  cex.fac = cex.fac)
    }
  }

  ## Also add a map
  make_map(obj1, tt, destin = destin)

  ## Also make a covariate plot
  plot_covariates(obj1, tt)

}





##' Plot a single cytogram.
##'
##' @param y (nt x 2) matrix.
##' @param counts multiplicity of each point in y.
##' @param cex_fac Only active when \code{!is.null(counts)}; user-supplier
##'   multiplier onto the point size \code{cex==sqrt(counts)}.
##'
##' @return NULL
one_2d_plot <- function(y, counts=NULL, xlim=NULL, ylim=NULL, xlab=NULL, ylab=NULL, cex=0.5,
                        pt_col = rgb(0, 0, 1, 0.1),
                        cex_fac = 1,
                        axes = TRUE,
                        x_ticks = NULL,
                        y_ticks = NULL){

  ## Basic checks.
  stopifnot(ncol(y) == 2)
  if(!is.null(counts)) stopifnot(length(counts) == nrow(y))

  if(is.null(xlim)) xlim = range(y[,1])
  if(is.null(ylim)) ylim = range(y[,2])

  ## Create empty plot
  plot(NA,
       ylim = ylim,
       xlim = xlim,
       ylab = ylab,
       xlab = xlab,
       cex.lab = 2,
       cex.axis = 2,
       xaxt = 'n',
       yaxt = 'n')
  if(!axes){
    axis(1, at = x_ticks,
         cex.axis = 2)
    axis(2, at = y_ticks,
         cex.axis = 2)
  } else {
    axis(1, cex.axis = 2)
    axis(2, cex.axis = 2)
  }

  ## Add datapoints
  if(is.null(counts)){
    cex = 0.01
  } else {
    cex = counts %>% sqrt()
    cex = cex * cex_fac
  }
  points(y, col = pt_col, pch = 16, cex = cex)
}



make_map_ggplot <- function(res, tt, destin=NULL){
  load(file = file.path("/home/shyun/Dropbox/research/usc/flow-cytometry/figures/chlorophyll-data.Rdata"),
       verbose = TRUE)
  ## sfl_small = sfl[c(500,1000,1500,2000,2500,3000,3500),]
  ## sfl_small

  ## ## Load lon/lat
  ## load(file = file.path("/home/shyun/repos/cruisedat/export/MGL1704-hourly-only-binned.Rdata"))
  ## X = X[-(1:12),]
  ## ylist = ylist[-(1:12)]
  ## countslist = countslist[-(1:12)]
  ## sfl_small = X[, c("lat", "lon")]
  ## save(sfl_small, file=file.path("/home/shyun/Dropbox/research/usc/flow-cytometry/figures/sfl_small.Rdata"))
  load(file.path("/home/shyun/Dropbox/research/usc/flow-cytometry/figures/sfl_small.Rdata"))

  ## Make the plot
  theme_set(theme_bw())
  world <- ne_countries(scale = "medium", returnclass = "sf")
  ggplot(data = world) +
    geom_point(data=dat, aes(x=lon, y=lat, color=SST), alpha=1, size = 20, shape=15, stroke=FALSE) +
    geom_point(data=dat, aes(x=lon, y=lat, color=SST), alpha=1, size = 7, shape=15, stroke=FALSE) +
    scale_colour_gradient(low = "darkblue", high = "red", limits = c(5, 28),
                          na.value = "red") +
    geom_sf() +
    xlab("") + ylab("") +
    coord_sf(xlim = xlim, ylim = ylim, expand = FALSE,
             label_axes = list(right = "N", bottom = "E")) +
    geom_point(dat=sfl_small[tt,], aes(x=lon, y=lat), size= 5) + ##, col='black')
    geom_line(dat=sfl_small, aes(x=lon, y=lat), size= 2) + ##, col='black')
    theme(axis.text = element_text(size=15),
          axis.title = element_text(size=15))
}
