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
