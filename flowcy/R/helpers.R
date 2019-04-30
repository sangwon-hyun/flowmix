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
                          fill=FALSE, beep=FALSE){

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
        if(beep & isim==nsim){beepr::beep()}
    }
    if(fill) cat(fill=TRUE)
}



##' Crude function to check convergence.
check_converge <- function(old, new, tol=1E-6){ return(abs(new-old) < tol)  }
check_converge_rel <- function(old, new, tol=1E-6){ return(abs((old-new)/old) < tol )  }


##' From applying k-means, ge tthe best cluster according to the gap statistic.
get_best_kmean_numclust <- function(data){
  ## data = iris[,1:3]
  gap_stat <- cluster::clusGap(data, FUN = kmeans, nstart = 25,
                      K.max = 10, B = 50, verbose=FALSE)##, method = "firstSEmax"
  with(gap_stat, cluster::maxSE(Tab[,"gap"],Tab[,"SE.sim"]))
}
