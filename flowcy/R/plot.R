##' Default printing functionality for covarem() results.
##' @param obj An object of class "covarem"
##' @return Nothing.
plot.covarem <- function(obj){

  ## Setup
  betas = do.call(cbind, lapply(obj$beta, function(beta)beta[-1,]))
  alpha =obj$alpha[ ,-1]
  cols = RColorBrewer::brewer.pal(obj$numclust, "Set3")
  dimdat = obj$dimdat

  ## Plot means, in each dimension
  par(ask=TRUE)
  for(idim in 1:dimdat){
    mn = obj$mn[,idim,]
    par(mfrow=c(1,2))
    matplot(mn, type='l', lty=1, lwd=2, main=paste0("Means in ", idim, "-th dim"),
            col=cols)
    wt = obj$pie
    matplot(wt, type='l', lty=1, lwd=2, main="Weight", ylim=c(0,1),
            col=cols)
  }

  ## Plot objective values
  plot(res$obj[-1], type='o', main="Objectives", pch=16, ylab="value", xlab="iterations")
  par(ask=FALSE)
}

plot2d <- function (x, ...) {
   UseMethod("plot2d", x)
}

##' Only in the 2d case, make series of plots.
##' @param show.fewer NULL by default. Otherwise the indices of the time points
##'   whose plots are to be shown.
##' @param obj Result of running covarem().
##' @param ylist List of responses.
##' @return Nothing
plot2d.covarem <- function(obj, ylist, ask=TRUE, show.fewer=NULL,
                           show.xb.constraint=FALSE
                           ## maxdev=NULL## Temporary addition.
                           ){

  ## Basic checks
  assert_that(obj$dimdat==2)

  ## Collect plot ranges
  all.y = do.call(rbind, ylist)
  ylim = range(all.y[,2])
  xlim = range(all.y[,1])

  ## ## Make plot
  ## par(mfrow=c(1,1))
  ## par(ask=TRUE)
  ## for(tt in 1:TT){
  ##   mn = obj$mn[tt,,]
  ##   wt = obj$pie[tt,]
  ##   plot(NA, ylim=ylim, xlim=xlim)
  ##   points(ylist[[tt]], col='grey50', pch=16)
  ##   points(mn[1,]~mn[2,],
  ##        cex=wt*10, pch=16, col='red')
  ## }
  ## Extract numclust  (eventually from somewhere else)
  mns = obj$mn
  numclust = obj$numclust
  TT =obj$TT
  p = ncol(obj$X)

  ## General plot settings
  cols = RColorBrewer::brewer.pal(numclust, "Set3")

  ## Eventually extract TT some other way.
  if(!is.null(show.fewer)){
    TTrange = show.fewer
    assert_that(all(TTrange <= TT))
    assert_that(all(TTrange >= 1))
  } else { TTrange = 1:TT }

  if(ask)par(ask=TRUE)
  par(mfrow=c(1,3))
  for(tt in TTrange){
    main0 = paste0("time ", tt, " out of ", TT)
    plot(NA, ylim=ylim, xlim=xlim, cex=3, ylab="", xlab="", main=main0)

    ## Add datapoints
    points(ylist[[tt]], col='grey90', pch=16, cex=.5)



    ## Begin temporary
    if(show.xb.constraint){

      ## (Temporary addition) Plot together all the centers.
      for(kk in 1:numclust){
        points(x=mns[1:TT,1,kk],
               y=mns[1:TT,2,kk], col='grey60', pch=16, cex=0.5)
      }

      ## Also add circle around beta0kwhose radius to maxdev.
      for(kk in 1:numclust){
        beta0list = lapply(obj$beta, function(betamat){
          betamat[1,]
        })
        points(x=beta0list[[kk]][1],
                y=beta0list[[kk]][2], col="blue", pch=16)
        print(obj$maxdev)
        plotrix::draw.circle(x=beta0list[[kk]][1],
                             y=beta0list[[kk]][2],
                             radius=obj$maxdev,
                             border="blue", lwd=2)
      }
    }
    ## End of temporary

    ## Collect pies
    pies = lapply(1:numclust, function(iclust){  obj$pie[,iclust] })

    ## Add fitted means
    for(kk in 1:numclust){
      points(x=mns[tt,1,kk],
             y=mns[tt,2,kk],
             col='red', pch=1, cex=pies[[kk]][tt]*10)

      text(x=mns[tt,1,kk],
           y=mns[tt,2,kk],
           labels=kk,
           col='red', pch=16, cex=pies[[kk]][tt]*5)
    }

    ## Add legend
    legend("bottomright",
           col=c("black", "red", "blue", "grey40"),
           pch=c(16,16,16, 16),
           legend = c("truth", "fitted", "intercept", "other fitted"))

    for(kk in 1:numclust){
    lines(ellipse::ellipse(x=obj$sigma[tt,kk,,],
                           centre=mns[tt,,kk]
                           ), lwd=1/2, col='red', lty=2)
    }
    ## lines(ellipse(rho), col="red")       # ellipse() from ellipse package
    ## lines(ellipse(rho, level = .99), col="green")
    ## lines(ellipse(rho, level = .90), col="blue")


    ## Plot pies
    plot(NA, xlim=c(0,TT), ylim=c(0,1.4), main=main0,
         ylab = "Mixture probability", xlab="time, t=1,..,T")
    for(iclust in 1:numclust){
      lines(pies[[iclust]], col=cols[iclust], lwd=2, type='o', pch=toString(iclust))
    }
    abline(v=tt, col='green')


    ## Plot covariates
    ylim.cov=range(obj$X)*1.5
    plot(NA, xlim=c(0,TT), ylim=ylim.cov, main="Covariates",
         ylab = "Environmental covariates", xlab="time, t=1,..,T")
    for(ii in p:1){
      lines(obj$X[,ii], col=ii, lwd=2, type='l')
    }
    abline(v=tt, col='green')

    ## if(plot.iter>1){
      ## betas = obj$beta.list[[plot.iter]]
      ## betas = do.call(cbind, lapply(obj$beta, function(beta)beta[-1,]))
      betas = round(do.call(cbind, obj$beta),2)
      rownames(betas) = paste("beta:", c("intrcpt", paste0("coef", 1:5)))
      colnames(betas)= paste0("clust-", c(1,1,2,2,3,3,4,4))
      text(0.6,max(ylim.cov)/1.5, paste(capture.output(betas), collapse='\n'), pos=4)##, family="monospace")

      alphas = as.matrix(t(obj$alpha))
      rownames(alphas) = paste("alpha:", c("intrcpt", paste0("coef", 1:5)))
      alphas = round(alphas,2)
      colnames(alphas)= paste0("cluster ", c(1:4))
      text(25,min(ylim)/1.5, paste(capture.output(alphas), collapse='\n'), pos=4)##, family="monospace")
  }
  par(ask=FALSE)
}


## ##' Visualizes the  data only. Temporarily useful.
## make_2d_plot_of_data <- function(figdir, plotname, res, plot.iter=NULL, truths=NULL, h=8){

##   if(is.null(plot.iter)) plot.iter = res$final.iter
##   ## pdf(file=file.path(figdir, plotname), width=15, height=8)
##   pdf(file=file.path(figdir, plotname), width=h, height=h)

##   ## Extract numclust  (eventually from somewhere else)
##   mns = res$mn.list[[plot.iter]]
##   numclust = res$numclust


##   ## Define range of date
##   xlim = range(do.call(rbind, res$ylist)[,1])
##   ylim = range(do.call(rbind, res$ylist)[,2])

##   ## General plot settings
##   cols = RColorBrewer::brewer.pal(numclust, "Set3")

##   ## Eventually extract TT some other way.
##   for(tt in 1:res$TT){
##     main0 = paste0("time ", tt, " out of ", TT)
##     main = paste0("Iteration ", plot.iter)
##     plot(NA, ylim=ylim, xlim=xlim, cex=3, ylab="", xlab="", main=main0)

##     ## Add datapoints
##     points(res$ylist[[tt]], col='lightgrey', pch=16, cex=.5)

##     ## Add truth as well
##     if(!is.null(truths)){
##       numclust.truth = dim(truths$mns)[3]
##       for(kk in 1:numclust.truth){
##         points(x=truths$mns[tt,1,kk],
##                y=truths$mns[tt,2,kk],
##                col="black", pch=16, cex=truths$pies[[kk]][tt]*5)
##       }
##     }

##   }
##   graphics.off()
## }
