##' Plot 1 dimensional model results.
##' @param ylist Data.
##' @param countslist Multiplicity.
##' @param res Optionally, |covarem| object of fitted model.
##' @param scale Defaults to TRUE. If TRUE, scale the data colors, relative to
##'   the largest count in |countslist|.
##' @param main Plot title
##'
##' @return NULL.
##' @export
plot_1d <- function(ylist, countslist, res = NULL, scale = TRUE, main = "",
                    date_axis = TRUE, no_axis = FALSE, mn.scale = 5,
                    cex_clust_label = 1.5,
                    omit_band = FALSE,
                    omit_label = FALSE,
                    reorder_clusters = TRUE,
                    cex_data=15){

  ## Plot ylist only
  stopifnot(ncol(ylist[[1]]) == 1)
  stopifnot(ncol(countslist[[1]]) == 1)
  plot_ylist(ylist, countslist, scale = scale, main = main, cex=cex_data)

  ## Plot model if applicable.
  if(!is.null(res)){
    stopifnot(res$dimdat==1)
    cols = RColorBrewer::brewer.pal(max(3,res$numclust), "Set2")
    if(res$numclust>8){
      cols = RColorBrewer::brewer.pal(max(3,res$numclust), "Paired")
    }

    ## Reorder clusters, according to total pie.
    if(reorder_clusters) res <- reorder_clust(res)

    ## Plot the means
    add_mn(res, cols, mn.scale)

    ## Add confidence band
    if(!omit_band){  add_band(res, cols)   }

    ## Make the remaining region outside of TT white, to prevent spillover of
    ## the cluster means' cex=17 points
    fill_both_sides(res$TT)

    ## Make the X ticks dates.
    if(!no_axis){
      if(date_axis){ add_date_ticks(res) } else { axis(1) }
    }

    ## Add Cluster labeling
    if(!omit_label){
    add_cluster_label(res, cex_clust_label)
    }
  }
  return(cols)
}

##' (incomplete)
plot_X <- function(res){

  ## 3. Plot the covariates
  matplot(X, type='l', lty=1, ylab="", xlab="Time", lwd=1, col="grey50", axes=FALSE)
  cols = RColorBrewer::brewer.pal(max(3,ncol(X)), "Set3")
  axis(1); axis(2);

  title(main="Covariates", cex.main=2)
  ## colnames(X)[1:3] = c("b1", "b2", "b2")## temporary
  if(!is.null(colnames(X)))legend("bottomright", lwd=1, col=cols, legend=colnames(X), cex=1.3,
                                  bg="white",  ncol=2)
}

##' Make 1d data plot.
plot_ylist <- function(ylist, countslist, res = NULL, scale = TRUE, main = "",
                       cex = 1){
  stopifnot(ncol(ylist[[1]]) == 1)
  TT = length(ylist)
  matplot(NA, type='l', lty=1, lwd=.1,
          ylim = range(unlist(ylist)),
          xlim = c(1, TT),
          ylab="", xlab="",
          axes=FALSE)

  ## Add title
  title(main = main)

  ## Handle when countslist is NULL
  if(is.null(countslist)) countslist = sapply(ylist, function(y)rep(1,nrow(y)))

  ## Visualize the data
  if(!scale) mx = 10 else mx = max(unlist(countslist))
  for(tt in 1:TT){
    y = ylist[[tt]]
    ct = countslist[[tt]] / mx
    points(x = rep(tt, length(y)),
           y = y, col = rgb(0, 0, 0, ct),
           pch = 15, cex = cex)
  }
}

##' Add cluster label near beginning of plot (only for 1d data).
##' @param cex if NULL, don't make cluster label.
add_cluster_label <- function(res, cex=1.5){
  if(!is.null(cex)){
    for(iclust in 1:res$numclust){
      y = res$mn[1, 1, iclust]
      text(x = res$TT/50, y = y,
           label = paste0("Cluster ", iclust), cex = cex)
    }
  }
}

##' Add mean (only for 1d data)
add_mn <- function(res, cols, mn.scale = 5){
  for(iclust in 1:res$numclust){
    lines(res$mn[,1,iclust], type = 'o', lty = 1, lwd = .1,
             cex = res$pie[,iclust] * mn.scale, pch = 15, col = cols[iclust])
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
    polygon(c(1:length(up),rev(1:length(up))),
            c(up,rev(dn)),col=grDevices::adjustcolor(cols[iclust], alpha.f = 0.2),
            border=NA)
  }
}

##' Make the remaining region outside of TT white, to prevent spillover of the
##' cluster means' cex=17 points. (only for 1d data).
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

##' Make ticks from rownames of res$X. TODO: make it handle dates. (only for 1d
##' data).
##'
##' @param res Object of class |covarem|. The row names of res$X are used.
##' @param ... Rest of arguments to both axes via \code{axis()}.
##'
##' @return No return.
add_date_ticks <- function(res, ...){
  add_date_ticks_from_dates(rownames(res$X), ...)
  ## dates = sapply(as.Date(rownames(res$X)) %>% format("%B %d"), toString) ## OR manually bring in,
  ##                                                      ## using an argument to
  ##                                                      ## the function.
  ## nums = as.numeric(as.factor(dates))
  ## left_ticks = sapply(sort(unique(nums)),function(ii){min(which(nums==ii))})
  ## left_ticks = c(left_ticks, res$TT)
  ## mid_ticks = sapply(sort(unique(nums)),function(ii){mean(which(nums==ii))})
  ## dates_mid_ticks = dates[round(mid_ticks)]
  ## axis(1, at=left_ticks, labels=FALSE)
  ## axis(1, at=mid_ticks, labels = dates_mid_ticks, tick=FALSE, las=2, ...)
  ## axis(2, ...)
}


##' Add date ticks from string of dates.
##'
##' @param dates Vector of strings of the form "2017-05-29T00:00:00".
##' @param ... Rest of arguments to both axes via \code{axis()}.
##'
##' @return No return.
add_date_ticks_from_dates <- function(dates, ...){
  dates = sapply(as.Date(dates) %>% format("%B %d"), toString) ## OR manually bring in,
                                                       ## using an argument to
                                                       ## the function.
  nums = as.numeric(as.factor(dates))
  left_ticks = sapply(sort(unique(nums)),function(ii){min(which(nums==ii))})
  left_ticks = c(left_ticks, length(dates))##res$TT)
  mid_ticks = sapply(sort(unique(nums)),function(ii){mean(which(nums==ii))})
  dates_mid_ticks = dates[round(mid_ticks)]
  axis(1, at=left_ticks, labels=FALSE)
  axis(1, at=mid_ticks, labels = dates_mid_ticks, tick=FALSE, las=2, ...)
  axis(2, ...)
}



##' Reorder the results so that cluster 1 through numclust is in decreasing
##' order of the total probability, over all time points.
reorder_clust <- function(res){

  ## Create an order

  ## Here's a suggestion for the 1d plots... since the coloring (and
  ## likewise numbering) of the clusters is arbitrary, what if we come up
  ## with some standard rule for labeling/coloring.  For example, it could
  ## simply be from largest diameter (averaged across all time) to smallest
  ## diameter.  Another natural choice for the ordering would be to order
  ## them from largest to smallest pi values (again averaged over all time).
  ## This way, your 1d-all-models.pdf will have (for the most part)
  ## consistent coloring/labeling in all the tiny plots.  Actually, the
  ## "largest pi" approach would work for the 3d plots as well.

  ## ord = res$pie %>% colSums %>% order(decreasing=TRUE)
  ## ord = res$mn %>% colSums %>% order(decreasing=TRUE)
  ord = res$mn[,1,] %>% colSums %>% order(decreasing=TRUE)

  ## Reorder mean
  res$mn = res$mn[,,ord, drop=FALSE]

  ## Reorder sigma
  res$sigma = res$sigma[ord,,,drop=FALSE]

  ## Reorder pie
  res$pie = res$pie[,ord, drop=FALSE]

  ## Reorder the alpha coefficients
  res$alpha = res$alpha[ord,, drop=FALSE] ## Also rename the row names
  rownames(res$alpha) = paste0("clust-", 1:res$numclust)

  ## Reorder the beta coefficients
  res$beta = res$beta[ord]
  names(res$beta) = paste0("clust-", 1:res$numclust)

  return(res)
}
