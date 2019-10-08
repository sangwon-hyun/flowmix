##' A helper function to print the progress of a simulation.
##' @param isim isim.
##' @param nsim number of sim.
##' @param type type of job you're running. Defaults to "simulation".
##' @param lapsetime lapsed time, in seconds (by default).
##' @param lapsetimeunit "second".
##' @param start.time start time.
##' @param fill Whether or not to fill the line.
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


##' In  an "opp"  table, remove  all the  top-coded data  points (rows)  in each
##' variable.
##' @param opp opp table.
##' @return Cleaned version of \code{opp} with topcoded rows removed.
##' @export
remove_topcoded <- function(opp){

  ## Remove how many are top-coded from below.
  which_chl_topcode = which(log(opp[,c("chl_small")])==0)
  which_pe_topcode = which(log(opp[,c("pe")])==0)
  which_fsc_topcode = which(log(opp[,c("fsc_small")])==0)

  ## Remove all top coded from above.
  which_chl_topcode2 = which((opp[,c("chl_small")])==max(opp[,c("chl_small")]))
  which_pe_topcode2 = which((opp[,c("pe")])==max(opp[,c("pe")]))
  which_fsc_topcode2 = which((opp[,c("fsc_small")])==max(opp[,c("fsc_small")]))

  opp = opp[-unique(c(which_chl_topcode,
                      which_pe_topcode,
                      which_fsc_topcode,
                      which_chl_topcode2,
                      which_pe_topcode2,
                      which_fsc_topcode2
                      )),]
  return(opp)
}

##' A temporary plotting function.
##' @param mytable opp table to plot.
##' @param chl_topcode The number of datapoints whose "CHL" is topcoded.
##' @param pe_topcode The number of datapoints whose "PE" is topcoded.
##' @param fsc_topcode The number of datapoints whose "FSC" is topcoded.
##' @param col Color of points.
##' @param ... Additional arguments to \code{plot()}.
myplot <- function(mytable, chl_topcode=NULL,
                   pe_topcode=NULL, fsc_topcode=NULL,
                   col = "#00000055",...){

  main = paste0(nrow(mytable),
                " robust particles, log scale \n",
                paste0("topcoded (at 1): chl=", chl_topcode, ", ",
                       "pe=", pe_topcode, ", ",
                       "fsc=", fsc_topcode))
  graphics::plot(mytable, col=col,
       main=main,
       pch=16,
       cex=1,
       ...)
}



## Helper function to logarithmically space out R.  |length| values linear on
## the log scale from |to| down to |from|.  TODO: Write this in:
## "lambda.min.ratio=ifelse(nobs<nvars,1e-2,1e-4)".
logspace <- function(max, min=NULL, min.ratio = 1E-4, length){
  if(is.null(min)) min = max * min.ratio
  matlab::logspace(log10(min), log10(max), length)
}
