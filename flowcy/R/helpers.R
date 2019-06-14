##' A helper function to print the progress of a simulation.
##' @param isim isim.
##' @param nsim number of sim.
##' @param type type of job you're running. Defaults to "simulation".
##' @param lapsetime lapsed time, in seconds (by default).
##' @param lapsetimeunit "second".
##' @param start.time start time.
##' @param fill Whether or not to fill the line.
##' @param beep Whether to beep when done.
printprogress <- function(isim, nsim, type="simulation", lapsetime=NULL,
                          lapsetimeunit="seconds", start.time=NULL,
                          fill=FALSE){

    ## If lapse time is present, then use it
    if(fill) cat(fill=TRUE)
    if(is.null(lapsetime) & is.null(start.time)){
            cat("\r", type, " ", isim, "out of", nsim)
    } else {
        if(!is.null(start.time)){
            lapsetime = round(difftime(Sys.time(), start.time,
                                       units = "secs"), 0)
            remainingtime = round(lapsetime * (nsim-isim)/isim,0)
            endtime = Sys.time() + remainingtime
        }
        cat("\r", type, " ", isim, "out of", nsim, "with lapsed time",
            lapsetime, lapsetimeunit, "and remaining time", remainingtime,
            lapsetimeunit, "and will finish at", strftime(endtime))
        ## if(beep & isim==nsim){beepr::beep()}
    }
    if(fill) cat(fill=TRUE)
}



##' Crude function to check convergence.
check_converge <- function(old, new, tol=1E-6){ return(abs(new-old) < tol)  }
check_converge_rel <- function(old, new, tol=1E-6){ return(abs((old-new)/old) < tol )  }


##' From applying k-means, get the best cluster according to the gap statistic.
get_best_kmean_numclust <- function(data){
  gap_stat <- cluster::clusGap(data, FUN = kmeans, nstart = 25,
                      K.max = 10, B = 50, verbose=FALSE)##, method = "firstSEmax"
  with(gap_stat, cluster::maxSE(Tab[,"gap"],Tab[,"SE.sim"]))
}



##' @param obj is output from driftem
##' @param plotdir directory to plot to "~/Desktop/cytometry-plots"
##' @param filedir name of plot
plot_driftem <- function(res, plotdir="~/Desktop", plotname="driftem-example.pdf",
                         width=10, height=10, dimplot=1:2){


  ## For now, only written for 2d images. If not 2d image, then nothing is plotted.
  if(obj$dimdat!=2){
    warning(paste0("only the first two dimensions ",
                   "(or the dimensions you prescribed in |dimplot| are being plotted!"))
  }
  if(length(dimplot)!=2) stop("only two dimensions can be plotted at a time")

  ## Setup
  alldat = do.call(rbind, obj$data)
  xlim = range(alldat[,dimplot[1]])
  ylim = range(alldat[,dimplot[2]])

  pdf(file=file.path(plotdir, plotname), width=width, height=height)

  ## For every time point
  for(tt in 1:(res$TT)){
    ## if(tt==2)browser()
    printprogress(tt, (res$TT), "plots")

    plot(data[[tt]][,dimplot], xlim=xlim, ylim=ylim)## , col=labels[[tt]]

    ## Plot fitted means
    points(res$mulist[[res$final.iter]][tt,,],
           cex=res$pielist[[res$final.iter]][tt,]*10,
           pch=16, col=c(1,2,3))

    ## Plot fitted covariance contours
    all.lims = cbind(ylim,xlim)##sapply(1:2, function(mycol){range(data[[tt]][,mycol])})
    for(jj in 1:numclust){
      add.contours(all.lims, c(1,2),
                   res$mulist[[res$final.iter]],
                   res$sigmalist[[res$final.iter]],
                   tt, jj, lwd=.5)
    }
  }
  graphics.off()
}


##' Obtain spectral range from |res|.
spectral_range_from_res <- function(res, fac=2){
  ## warmstart.muarray = init_mu_warmstart_v2(data, numclust, TT, 1E-6, verbose=FALSE)
  alldata = do.call(rbind, res$data)
  res = driftem(alldata,
                mu=res$init_mu, ##warmstart.muarray,
                numclust=res$numclust)
  sigma = res$sigmalist[[res$final.iter]][1,,,]
  sigmalist = lapply(1:dim(sigma)[1], function(idim){sigma[idim,,]})
  mydecomp = function(sigma, fun){
    decomp = eigen(sigma)
    eigvals = decomp$values
    return(fun(eigvals))
  }
  return(list(min = min(sapply(sigmalist, mydecomp, min)) * (1/fac),
              max = max(sapply(sigmalist, mydecomp, max)) * fac))
}
