##' Let's make an artificial dataset of 3 images, from iris
make_data <- function(T, fac=1, noise = 0.2){
  set.seed(0)
  data(iris)
  data = list()
  data[[1]] = iris[,1:3]
  ## data[[2]] = iris[,1:3] + fac*rbind(c(0.5, -0.5, 0.5)) + mvrnorm(n=nrow(iris), mu=rep(0,3),
  ##                                                          Sigma = diag(rep(0.2,3)))
  ## data[[3]] = iris[,1:3] + fac*rbind(c(0.5, -0.5, 0.5)) + mvrnorm(n=nrow(iris), mu=rep(0,3),
                                                           ## Sigma = diag(rep(0.2,3)))

  data[[2]] = iris[,1:3] + fac*rbind(c(1,1,1)) + mvrnorm(n=nrow(iris), mu=rep(0,3),
                                                           Sigma = diag(rep(noise,3)))
  data[[3]] = iris[,1:3] + 2*fac*rbind(c(1,1,1)) + mvrnorm(n=nrow(iris), mu=rep(0,3),
                                                           Sigma = diag(rep(noise,3)))
  stopifnot(T==3)
  return(data[1:T])
}



##' Helper
myplot <- function(data, mu, pie=NULL, sigma, lam1=NULL, lam2=NULL, T=3, dimdat=3, numclust=2, fac=8, lwd=1, cex.data=1){

  stopifnot(dimdat==3)

  ## Extract the axis limits
  iis.list = list(c(1:2), (2:3), c(3,1))
  all.dims = 1:3
  all.lims = lapply(all.dims, function(mydim){
    sapply(1:3, function(mytime){
      (data[[mytime]][,mydim])
    })
  })
  all.lims = do.call(rbind, all.lims)
  all.lims = apply(all.lims, 2, range)

  par(mfrow=c(3,3))
  for(t in 1:T){
    for(iis in iis.list){
      plot(data[[t]][,iis], pch=16, main=paste0("time=",t, " lam1=", lam1, ", lam2=", lam2),
           xlim = all.lims[,iis[2]], ylim = all.lims[,iis[1]], cex=cex.data)
           ##, xlim=c(0,10), ylim=c(0,10))

      ## Add cluster centers
      if(!is.null(pie)) cex = pie[t,]*fac else cex=2
      points(mu[t,,iis], pch=16, col='red', cex=cex)

      ## Add *other* cluster centers
      for(t.other in (1:T)[-t]){
        if(!is.null(pie)) cex = pie[t.other,]*8 else cex=2
        points(mu[t.other,,iis], pch=16, col=rgb(0,1,0,0.5), cex=cex)
      }

      ## Add contours
      for(jj in 1:numclust){
        add.contours(all.lims, iis, mu, sigma, t, jj, lwd=lwd)
      }
    }
  }
}


##' Helpers to contours, to myplot().
add.contours <- function(all.lims, iis, mu, sigma, t, jj, lwd=1, col='grey'){
  x.points = seq(from=all.lims[1,iis[2]], to=all.lims[2,iis[2]], length=100)
  y.points = seq(from=all.lims[1,iis[1]], to=all.lims[2,iis[1]], length=100)
  z <- matrix(0,nrow=100,ncol=100)
  for (i in 1:100) {
    for (j in 1:100) {
      z[i,j] <- mvtnorm::dmvnorm(c(x.points[i],y.points[j]),
                        mean=mu[t,jj,iis],sigma=sigma[t,jj,iis,iis])
    }
  }
  contour(x.points,y.points,z, add=TRUE,col=col, lwd=lwd)
}



##' Computes x = A^+ * b using an SVD (slow but stable)
svdsolve <- function(A,b,rtol=1E-20) {
  s = svd(A)
  di = s$d
  ii = di>rtol
  di[ii] = 1/di[ii]
  di[!ii] = 0
  return(list(x=s$v%*%(di*(t(s$u)%*%b)),q=sum(ii)))
}



##' plotter for *series* of plots.
myplots <- function(plotname=NULL, data, labels, mns1, mns2, mns3, sigmas1,
                    sigmas2, sigmas3, results){
  if(is.null(plotname))plotname = "example-dat.pdf"
  pdf(file=file.path("~/Desktop/cytometry-plots", plotname),
      width=10, height=10)
file.path("~/Desktop/cytometry-plots", plotname)
  for(tt in 1:TT){
    plot(data[[tt]][,1:2],
         col=labels[[tt]],
         xlim = c(-3,5), ylim=c(-3,6))

    ## Plot actual means
    for(ttt in (1:TT)[-tt]){
      points(rbind(mns1[[ttt]], mns2[[ttt]]), col="lightgrey",
             cex=2, pch="+")
    }
    points(rbind(mns1[[tt]], mns2[[tt]]), col=c(1,2), cex=3, pch="+")


    ## Plot fitted means
    points(results$mulist[[results$final.iter]][tt,,], col=c(1,2,3),
           cex=results$pielist[[results$final.iter]][tt,]*10,
           pch=16)

    ## Plot the covariances as well.
  }
  graphics.off()
}
