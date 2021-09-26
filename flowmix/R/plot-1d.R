##' Plot 1 dimensional model results.
##'
##' @param ylist Data.
##' @param countslist Multiplicity.
##' @param res Optionally, |covarem| object of fitted model.
##' @param scale Defaults to TRUE. If TRUE, scale the data colors, relative to
##'   the largest count in |countslist|.
##' @param main Plot title.
##' @param time Vector of date strings.
##'
##' @return NULL.
##' @export
plot_1d <- function(ylist,
                    countslist,
                    res = NULL,
                    scale = TRUE,
                    main = "",
                    time_axis = FALSE,
                    time = NULL,
                    mn.scale = 5,
                    cex_clust_label = 1.5,
                    omit_band = FALSE,
                    omit_label = FALSE,
                    reorder_clusters = TRUE,
                    cex_data = 1,##15
                    ylim = NULL,
                    cols = NULL){

  ## Plot ylist only
  stopifnot(ncol(ylist[[1]]) == 1)
  stopifnot(ncol(countslist[[1]]) == 1)
  plot_ylist(ylist, countslist, scale = scale, main = main, cex=cex_data, ylim = ylim)
  if(!is.null(res)) assertthat::assert_that(is(res, "flowmix"))

  ## Make the X ticks dates.
  if(time_axis){

    ## Find time.
    if(!is.null(time)){
      assertthat::assert_that(check_if_date(time))
    } else {
      if(!is.null(res)){
        check_if_date(rownames(res$X))
        assertthat::assert_that(check_if_date(rownames(res$X)))
        time = res$X %>% rownames() %>% lubridate::as_datetime()
      }
    }
    add_date_ticks_from_dates(time)
  } else {
    axis(1)
    axis(2)
  }

  ## Plot model means, if provided.
  if(!is.null(res)){
    stopifnot(res$dimdat==1)

    if(is.null(cols)){
    if(res$numclust <= 8){
      cols = RColorBrewer::brewer.pal(max(3,res$numclust), "Set2")
    } else {
      cols = RColorBrewer::brewer.pal(max(3,res$numclust), "Paired")
    }
    }

    ## Reorder clusters, according to total prob.
    if(reorder_clusters) res <- reorder_clust(res)

    ## Plot the means
    add_mn(res, cols, mn.scale)

    ## Add confidence band
    if(!omit_band){  add_band(res, cols)   }

    ## Make the remaining region outside of TT white, to prevent spillover of
    ## the cluster means' cex=17 points
    fill_both_sides(length(ylist))

    ## Add Cluster labeling
    if(!omit_label){
      add_cluster_label(res, cex_clust_label)
    }
  }
  return(list(res = res, cols = cols))
}

##' Make 1d data plot. No axes are added.
##'
##' @inheritParams plot_1d
##'
##' @export
plot_ylist <- function(ylist, countslist, res = NULL, scale = TRUE, main = "",
                       cex = 1, ylim = NULL){

  stopifnot(ncol(ylist[[1]]) == 1)
  TT = length(ylist)
  if(is.null(ylim)) ylim = range(unlist(ylist))
  matplot(NA,
          type = 'l',
          lty = 1,
          lwd = .1,
          ylim = ylim,
          xlim = c(1, TT),
          ylab = "",
          xlab = "",
          axes = FALSE)

  ## Add title
  title(main = main)

  ## Handle when countslist is NULL
  if(!is.null(countslist)){
    if(!scale) mx = 10 else mx = max(unlist(countslist))
    cols = lapply(countslist, function(counts){
      ct = counts / mx
      return(rgb(0, 0, 0, ct))
    })
    pch = 15
  } else {
    countslist = sapply(ylist, function(y) rep(1, nrow(y)))
    cols = lapply(1:TT, function(a) rgb(0,0,0,0.1))
    pch = 16
  }

  ## Visualize the data
  for(tt in 1:TT){
    y = ylist[[tt]]
    col = cols[[tt]]
    points(x = rep(tt, length(y)),
           y = y, col = col,
           pch = pch, cex = cex)
  }
}


##' Add cluster label near beginning of plot (only for 1d data).
##' @param cex if NULL, don't make cluster label.
##' @export
add_cluster_label <- function(res, cex=1.5){
  if(!is.null(cex)){
    for(iclust in 1:res$numclust){
      y = res$mn[1, 1, iclust]
      graphics::text(x = res$TT/50, y = y,
                     label = paste0("Cluster ", iclust), cex = cex)
    }
  }
}

##' Add mean (only for 1d data)
add_mn <- function(res, cols, mn.scale = 5){
  for(iclust in 1:res$numclust){
    graphics::lines(res$mn[,1,iclust], type = 'o', lty = 1, lwd = .1,
                    cex = res$prob[,iclust] * mn.scale, pch = 15, col = cols[iclust])
  }
}

##' Add bands (only for 1d data).
add_band <- function(res, cols){
  sigmas = sqrt(res$sigma[,1,])
  sigma_mat = matrix(sigmas,## res$sigma[,1,],
                     ncol=res$numclust,
                     nrow=res$TT, byrow=TRUE)
  for(iclust in 1:res$numclust){
    up = res$mn[,1,iclust] + 2*sigma_mat[,iclust]
    dn = res$mn[,1,iclust] - 2*sigma_mat[,iclust]
    graphics::polygon(c(1:length(up),
                        rev(1:length(up))),
                      c(up,rev(dn)),
                      col = grDevices::adjustcolor(cols[iclust], alpha.f = 0.2),
                      border = NA)
  }
}

##' Make the remaining region outside of TT white, to prevent spillover of the
##' cluster means' cex=17 points. (only for 1d data).
##'
fill_both_sides <- function(TT){
  x = c(TT -1 + (1:100), TT -1 + (100:1))
  polygon(x = x, y = c(rep(-5, 100), rep(5, 100)),
          col = "white",
          border=FALSE)

  x = c(3-(1:100), 3-(100:1))
  polygon(x = x, y = c(rep(-5, 100), rep(5, 100)),
          col = "white",
          border=FALSE)
}


##' Check if time is TRUE.
##'
##' @param date String vector containing date strings
##'   e.g. "2017-06-13T11:00:00".
##'
##' @return TRUE only all are dates..
check_if_date <- function(date){
  date = date %>% lubridate::as_datetime()
  if(length(date) == 0) return(FALSE)
  return(all(sapply(date, lubridate::is.POSIXt)))
}



##' Add date ticks from string of dates.
##'
##' @param dates Vector of strings of the form "2017-05-29T00:00:00", or
##'   otherwise recognizeable using \code{lubridate::as_datetime()}.
##' @param empty_tick_marks Defaults to FALSE. If TRUE, use EMPTY tick marks.
##' @param ... Rest of arguments to both axes via \code{axis()}.
##'
##' @return No return.
##' @export
add_date_ticks_from_dates <- function(dates, empty_tick_marks=FALSE, ...){
  ## dates = sapply(as.Date(dates) %>% format("%B %d"), toString)

  ## Get dates and X coordinates
  dates = sapply(lubridate::as_datetime(dates) %>% format("%B %d"), toString)
  nums = as.numeric(as.factor(dates))

  ## Form the tick locations.
  left_ticks = sapply(sort(unique(nums)),function(ii){min(which(nums==ii))})
  left_ticks = c(left_ticks, length(dates))##res$TT)
  mid_ticks = sapply(sort(unique(nums)),function(ii){mean(which(nums==ii))})
  dates_mid_ticks = dates[round(mid_ticks)]
  if(empty_tick_marks) dates_mid_ticks = rep("", length(dates_mid_ticks))

  ## Place those ticks
  graphics::axis(1, at = left_ticks, labels = FALSE)
  graphics::axis(1, at = mid_ticks, labels = dates_mid_ticks, tick = FALSE, las=2, ...)
  graphics::axis(2, ...)
}


##' Make a 1d binned data plot, using ggplot.
##'
##' @param ylist Must be numeric vectors containing bin centers.
##' @param qclist Corresponding QC in those bins.
##'
##' @return ggplot
bin_plot_1d <- function(ylist, qclist, col=NULL){

  if(is.null(col)) col = c("white", "black", "yellow", "red")

  stopifnot(all(sapply(ylist, length) == sapply(qclist, length)))

  combined_dat = lapply(names(ylist), function(onedatetime){
    tibble(y = ylist[[onedatetime]],
           qc = qclist[[onedatetime]],
           time = onedatetime)
  }) %>% bind_rows()
  combined_dat = combined_dat %>% mutate(time = lubridate::as_datetime(time))

  p = combined_dat %>% ggplot() + geom_raster(aes(x=time, y=y, fill=qc)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    scale_fill_gradientn(colours = col)
  return(p)
}
