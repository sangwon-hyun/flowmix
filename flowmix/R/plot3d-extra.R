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




##' Make map for the MGL1704 cruise.
make_map <- function(res, tt, lon, lat){##, destin=NULL){

  ## Make empty map
  lat1 = 10
  lat2 = 60
  lon1 = -180
  lon2 = -120
  margin = 15
  xrange = c(lon1, lon2) + margin * c(-1,1)
  yrange = c(lat1, lat2) + margin * c(-1,1)
  maps::map(database = 'world',
            ## xlim = c(-170, -110),
            ## ylim = c(10, 50),
            xlim = xrange,
            ylim = yrange,
            fill = T,
            col = 'grey',
            resolution = 0,
            bg = 'white',
            mar = c(1,1,2,1))

  ## Make
  if(!is.null(lat) & !is.null(lon)){
    ## load(file=file.path(destin, "latlontime.Rdata"))
    colfunc <- colorRampPalette(c("red", "blue"))
    ## cols = colfunc(length(lat))
    cols = "black"
    lines(y = lat, x = lon, pch = 16, cex = 0.1, col = cols, lwd = .5)
  }
## Add the point for time tt
  points(lon[tt], lat[tt], pch = "X", cex = 4, col = 'red')
}


##' Make ticks from rownames of res$X. TODO: make it handle dates. (only for 1d
##' data).
##' @param res Object of class |flowmix|.
add_date_ticks <- function(res){
  dates = sapply(as.Date(rownames(res$X)), ## %>% format("%B %d")
                 toString)
  ## OR manually bring in,
  ## using an argument to
  ## the function.

  nums = as.numeric(as.factor(dates))
  left_ticks = sapply(sort(unique(nums)),function(ii){min(which(nums==ii))})
  left_ticks = c(left_ticks, res$TT)
  mid_ticks = sapply(sort(unique(nums)),function(ii){mean(which(nums==ii))})
  dates_mid_ticks = dates[round(mid_ticks)]
  axis(1, at=left_ticks, labels=FALSE)
  axis(1, at=mid_ticks, labels = dates_mid_ticks, tick=FALSE, las=2)
  axis(2)
}
