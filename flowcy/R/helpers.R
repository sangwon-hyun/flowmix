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
