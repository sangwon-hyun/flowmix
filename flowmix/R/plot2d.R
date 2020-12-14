##' CAN'T handle the binned case (with \code{countslist}) yet.
##'
##' @importFrom magrittr %>%
##'
##' @export
##'
plot2d <- function(ylist, res, tt, drawslist = NULL, resp = NULL, iclust = NULL, cols = NULL){

  ## Setup
  stopifnot(ylist %>% .[[1]] %>% ncol() == 2)
  rngs = do.call(rbind, ylist) %>% apply(., 2, range)
  numclust = res$numclust
  if(is.null(iclust)){ numclusts = 1:numclust } else {numclusts = iclust}
  if(is.null(cols)) cols = RColorBrewer::brewer.pal(numclust, "Set3")
  if(!is.null(cols)) assertthat::assert_that(length(cols) == numclust)


  cartoon_points <- function(y, x, cex, col, pch = 16, ...){
    points(y=y, x=x, cex=cex*10+2, pch = pch, col="black",...)
    points(y=y, x=x, cex=cex*10+1, pch = pch, col = col)
  }


  ## Draw membership
  if(is.null(drawslist)){
    if(is.null(resp)) resp <- Estep(res$mn, res$sigma, res$prob, ylist = ylist,
                                     numclust = res$numclust, first_iter = TRUE)
    drawslist = draw_membership(resp)
  }
  draws = drawslist[[tt]]

  ## Make plot
  plot(NA, ylim = rngs[,2], xlim = rngs[,1], xlab="", ylab="")

  ## Add the data points
  for(iclust in numclusts) points(ylist[[tt]][which(draws[,iclust]==1),1:2],
                                  col =cols[iclust] %>% adjustcolor(alpha.f = 0.5),
                                  pch = 16)
  ## Add the cluster centers
  for(iclust in numclusts){
    cartoon_points(y = res$mn[tt,2,iclust],
                   x = res$mn[tt,1,iclust],
                   cex = res$prob[tt,iclust],
                   col = cols[iclust] %>% adjustcolor(alpha.f = 1))
    lines(ellipse::ellipse(x = res$sigma[iclust,,], centre = res$mn[tt,, iclust]),
          lwd = 2, col = cols[iclust], lty = 1)
  }
}
