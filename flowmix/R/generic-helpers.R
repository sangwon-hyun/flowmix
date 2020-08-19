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



## Functions to check convergence.
check_converge <- function(old, new, tol=1E-6){ return(abs(new-old) < tol)  }
check_converge_rel <- function(old, new, tol=1E-6){ return(abs((old-new)/old) < tol )  }
check_converge_rel_print <- function(old, new, tol=1E-6){print(abs((old-new)/old)); return(abs((old-new)/old) < tol )  }


##' In an "opp" table (see the \code{popcycle} R package for more details),
##' remove all the top-coded data points (rows) in each variable.
##'
##' @param opp opp table.
##'
##' @return Cleaned version of \code{opp} with topcoded rows removed.
##'
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


##' Helper function to logarithmically space out R.  |length| values linear on
##' the log scale from |to| down to |from|.  TODO: Write this in:
##' "lambda.min.ratio=ifelse(nobs<nvars,1e-2,1e-4)".
logspace <- function(max, min=NULL, length, min.ratio = 1E-4){
  if(is.null(min)) min = max * min.ratio
  return(10^seq(log10(min), log10(max), length = length))
}
