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

  ## Extract numclust  (eventually from somewhere else)
  mns = obj$mn
  numclust = obj$numclust
  TT = obj$TT
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
  ## par(mfrow=c(1,3))
  m = matrix(c(1, 1, 2, 2, 4, 4,
               1, 1, 3, 3, 5, 5),
             nrow = 2, ncol = 6, byrow=TRUE)
  layout(m)
  par(mar=c(3,2,4,1))
  for(tt in TTrange){


    ## 1. Add main plot of data and fitted means

    ## Create empty plot
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
      lines(ellipse::ellipse(x = obj$sigma[tt, kk,,],
                             centre = mns[tt,, kk]
                             ), lwd=1/2, col='red', lty=2)
    }

    ## 2. Pies right now.
    pies.right.now = sapply(1:numclust, function(iclust){pies[[iclust]][[tt]]})
    names(pies.right.now) = paste("Clust", 1:numclust)
    barplot(pies.right.now,
            col = cols,
            lwd = 1,
            ylim = c(0, 1),
            main = paste0("Cluster prob. at time ", tt)
            )##, lwd=2, type='o', pch=toString(iclust))

    ## 3. Plot pies
    plot(NA, xlim=c(0,TT), ylim=c(0,1.4), main="Cluster probs.",
         ylab = "Mixture probability", xlab="time, t=1,..,T")
    for(iclust in 1:numclust){
      lines(pies[[iclust]], col=cols[iclust], lwd=2, type='o', pch=toString(iclust))
    }
    abline(v=tt, col='green')


    ## 4. Add table
    add_table(obj)

    ## 5. Plot covariates
    ylim.cov=range(obj$X)*1.5
    plot(NA, xlim=c(0,TT*1.3), ylim=ylim.cov, main="Covariates",
         ylab = "Environmental covariates", xlab="time, t=1,..,T")
    lwd = c(3,3,0.5,0.5,0.5)
    for(ii in p:1){
      lines(obj$X[,ii], col=ii, lwd=lwd[ii], type='l')
    }
    abline(v=tt, col='green')
    legend("topright", col=1:p, lwd=lwd, lty=1, legend=c("Salinity",
                                                         "SST", "Iron",
                                                         "Phosphorus",
                                                         "Chlorophyll"))
  }
  par(ask=FALSE)
}

add_table <- function(obj){

  ylim.cov=range(obj$X)*1.5

  ## Create the problem
  plot(NA, xlim=c(0,100), ylim=c(0,100),
       axes=FALSE)##c(0,1.4))##
       ## ylab = "Mixture probability", xlab="time, t=1,..,T")

  ## Add table of betas
  betas = round(do.call(cbind, obj$beta),2)
  rownames(betas) = paste("beta:", c("intrcpt", paste0("coef", 1:5)))
  betas[,seq(from = 1, to = ncol(betas), by=2)]

  betas = do.call(cbind, lapply(obj$beta, function(mybeta){
    vec = rep(NA, 2*nrow(mybeta))
    vec[seq(from = 1, to = nrow(mybeta)*2, by=2)] = mybeta[,1]
    vec[seq(from = 2, to = nrow(mybeta)*2, by=2)] = mybeta[,2]
    return(vec)
  }))
  rownames(betas) = c("Intercept", "",
                      "Salinity", "",
                     "SST", "",
                     "Iron", "",
                     "Phosphorus", "",
                     "Chlorophyll", "")
  betas = signif(betas, 2)
  betas = cbind(rep(c(1,2),6), betas)
  colnames(betas) = c("Dim", paste0("Clust ", c(1,2,3,4)))

  ## colnames(betas)= paste0("clust-", c(1,1,2,2,3,3,4,4))
  ## text(1,1,
  ##      paste(capture.output(betas), collapse='\n'),
  ##      pos=4)##, family="monospace")
  abg <- matrix(c(rep("grey90",10),
                  rep("grey80",10)), nrow=12, ncol=5,
                byrow=TRUE)
  plotrix::addtable2plot(-10,0,betas,
                         display.rownames = TRUE,
                         ## vlines=TRUE,
                         ## hlines=TRUE,
                         bg=abg,
                         title="Coef. for cluster centers")

  ## Add table of alphas
  alphas = as.matrix(t(obj$alpha))
  alphas = round(alphas,2)
  abg <- matrix(c(rep("grey90",4),
                  rep("grey80",4)), nrow=6, ncol=4,
                byrow=TRUE)
  rownames(alphas) = c("Intercept",
                      "Salinity",
                     "SST",
                     "Iron",
                     "Phosphorus",
                     "Chlorophyll")
  colnames(alphas) =  paste0("Clust ", c(1,2,3,4))
  plotrix::addtable2plot(50,0,alphas,
                         display.rownames = TRUE,
                         ## vlines=TRUE,
                         ## hlines=TRUE
                         bg=abg,
                         title="Coef. for cluster probs")

  ## rownames(alphas) = paste("alpha:", c("intrcpt", paste0("coef", 1:5)))
  ## colnames(alphas)= paste0("cluster ", c(1:4))
  ## text(25,min(ylim)/1.5, paste(capture.output(alphas), collapse='\n'), pos=4)##, family="monospace")
}


##' Makes a 'fancy' plotly plot from algorithm output.
##' @param res result of running covarem().
##' @return plotly object.
fancyplot <- function(res, saveplot=FALSE, filename=NULL,
                      title="3D Scatter plot"
                      ){
  data = lapply(1:res$numclust, function(iclust){
    x = res$mn[,1,iclust]
    y = res$mn[,2,iclust]
    z = 1:(res$TT)
    c = rep(iclust, res$TT)
    s = res$pie[,iclust]*50
    dat = data.frame(x, y, z, c, s)
    names(dat) = paste0( c("x", "y", "z", "c", "s"), iclust)
    dat
  })
  data = do.call(cbind, data)

  ## Setup
  scene = list(camera = list(eye = list(x = -1.25, y = 1.25, z = 1.25)),
               zaxis = list(title = "Time"),
               xaxis = list(title = "Dim 1"),
               yaxis = list(title = "Dim 2"))

  ## Make plot device
  p <- plotly::plot_ly(data, x = ~x1, y = ~y1, z = ~z1, type = 'scatter3d',
                       mode = 'lines+markers',
                       line = list(width = 2, fill = ~c1, colorscale = 'Viridis'),
                       marker = list(size = ~s1, fill = ~c1,
                                     colorscale = 'Viridis'),
                       name="Cluster 1")  %>%
  plotly::add_trace(x = ~x2, y = ~y2, z = ~z2,
            line = list(width = 2, fill = ~c2, colorscale = 'Viridis'),
            marker = list(size=~s2, fill=~c2),
            name="Cluster 2")%>% ## ,
  plotly::add_trace(x = ~x3, y = ~y3, z = ~z3,
            line = list(width = 2, fill = ~c3, colorscale = 'Viridis'),
            marker = list(size=~s3, fill=~c3),
            name="Cluster 3")%>% ## ,
  plotly::add_trace(x = ~x4, y = ~y4, z = ~z4,
            line = list(width = 2, fill = ~c4, colorscale = 'Viridis'),
            marker = list(size=~s4, fill=~c4),
            name="Cluster 4") %>%
  plotly::layout(title = title,
         scene = scene)
         ## autosize = F, width = 500, height = 500)


  if(saveplot) htmlwidgets::saveWidget(as.widget(p), filename)
  return(p)
}



fancytable <- function(res, type=c("alpha", "beta")){
  tables = get_table(res)
  if(type=="alpha")return(plotmat(tables$betas))
  if(type=="beta")return(plotmat(tables$alphas))

}

plotmat <- function(mat){
  p <- plot_ly(
    type = 'table',
    header = list(
        values = c("Dim", colnames(mat)),
        line = list(width = 1, color = 'black'),
        fill = list(color = 'rgb(235, 100, 230)'),
        font = list(family = "Arial", size = 14, color = "white")
    ),
    cells = list(
      values = rbind(
        rownames(mat),
        t(as.matrix(unname(mat)))
      ),
      align = c('left', rep('center', ncol(mat))),
      line = list(color = "black", width = 1),
      fill = list(color = c('rgb(235, 193, 238)', 'rgba(228, 222, 249, 0.65)')),
      font = list(family = "Arial", size = 12, color = c("black"))
    ))
  return(p)
}

get_table <- function(res){
  ## Add table of betas
  betas = round(do.call(cbind, res$beta),2)
  rownames(betas) = paste("beta:", c("intrcpt", paste0("coef", 1:5)))
  betas[,seq(from = 1, to = ncol(betas), by=2)]
  betas = do.call(cbind, lapply(res$beta, function(mybeta){
    vec = rep(NA, 2*nrow(mybeta))
    vec[seq(from = 1, to = nrow(mybeta)*2, by=2)] = mybeta[,1]
    vec[seq(from = 2, to = nrow(mybeta)*2, by=2)] = mybeta[,2]
    return(vec)
  }))
  rownames(betas) = c("Intercept", "",
                      "Salinity", "",
                     "SST", "",
                     "Iron", "",
                     "Phosphorus", "",
                     "Chlorophyll", "")
  betas = signif(betas, 2)
  betas = cbind(rep(c(1,2),6), betas)
  colnames(betas) = c("Dim", paste0("Clust ", c(1,2,3,4)))

  ## Add table of alphas
  alphas = as.matrix(t(res$alpha))
  alphas = round(alphas,2)
  abg <- matrix(c(rep("grey90",4),
                  rep("grey80",4)), nrow=6, ncol=4,
                byrow=TRUE)
  rownames(alphas) = c("Intercept",
                      "Salinity",
                     "SST",
                     "Iron",
                     "Phosphorus",
                     "Chlorophyll")
  colnames(alphas) =  paste0("Clust ", c(1,2,3,4))
  return(list(betas=betas, alphas=alphas))
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
