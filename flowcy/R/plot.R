#######################
## Function generics ##
#######################

plot2d <- function (x, ...) {
   UseMethod("plot2d", x)
}
plot2d_generic <- function (x, ...) {
   UseMethod("plot2d_generic", x)
}


##########################
## Function definitions ##
##########################

plot_covariates <- function(obj, tt = NULL){
  X = obj$X

  par(mar = c(5,5,2,2),
      cex.axis=2, cex.lab=2)
  ylim.cov = range(X) * c(1, 0.7)

  Xsmall = X %>% as_tibble %>% select(sss, sst, par)
  cols = 1:5
  TT = nrow(X)
  matplot(X, col='grey', lwd=0.5, type='l', xlim=c(0,TT*1), ylim = ylim.cov, axes=FALSE,
          ylab="", xlab="")
  legend("topleft", legend="Environmental Covariates", cex=3,
         bty = "n")
  add_date_ticks(obj)
  matlines(Xsmall, col=cols, lwd = 3)
  if(!is.null(tt)) abline(v = tt, col='green', lwd=2)

  Xsmall_names = c("Salinity", "Temperature", "Sunlight")
  ## ttt = 290 ## temporarily, this is a fixed place to put the labels.
  ttt = TT * 0.98 %>% round
  text(x=ttt, y=Xsmall[ttt,], label=Xsmall_names,## colnames(Xsmall)
       cex=2)
}

##' "Prettifies" covarem object. Not really a plotting function.
prettify <- function(res, signif_digit=2){

  ## Reorder clusters in decreasing order of diam
  res <- reorder_clust(res)

  ## Prettify betas
  betamat = do.call(cbind, res$beta)
  all.zero.rows = which(apply(betamat, 1, function(myrow)all(abs(myrow)<1E-8)))
  if(length(all.zero.rows) > 0){
    betamat = betamat[-all.zero.rows,, drop=FALSE]
  }
  betamat = betamat %>% Matrix::Matrix(sparse=TRUE) %>% signif(signif_digit)
  colnames(betamat) = unlist(Map(function(a,b) paste0(a, ", ", b),
                                 paste0("clust-", rep(1:res$numclust, each=res$dimdat)),
                                 colnames(betamat)))

  ## Prettify alphas
  alphamat = res$alpha %>% t %>% Matrix::Matrix(sparse=TRUE) %>% signif(signif_digit)
  all.zero.rows = which(apply(alphamat, 1, function(myrow)all(abs(myrow)<1E-8)))
  if(length(all.zero.rows) > 0){
    alphamat = alphamat[-all.zero.rows,, drop=FALSE]
  }

  ## Return
  return(list(alphamat = alphamat,
              betamat = betamat))
}


##' Plot pies.
##' @param res covarem object.
##' @param iclusts Optionally, provide the cluster numbers to plot just those.
##'
##' @return NULL
plot_pie <- function(res, iclusts=NULL, main=NULL,
                     cols = NULL
                     ){

  ## Setup
  if(is.null(iclusts)) iclusts = c(1:res$numclust)

  ## Reorder clusters in decreasing order of diam
  res <- reorder_clust(res)

  if(is.null(cols)){
    cols = RColorBrewer::brewer.pal(res$numclust, "Set3")
  }

  matplot(NA, xlim = c(0, res$TT),## res$pie
          ylab="", xlab="", ylim = c(0,1), axes = FALSE)
  abline(h = seq(from = 0, to = 1, by = 0.1), col='grey90', lwd=2, lty=3)
  matlines(res$pie[,iclusts], type='l', lty=1, lwd=3, col=cols[iclusts])
  ## axis(1); axis(2)
  if(is.null(main)){
    title(main="Cluster probabilities", cex.main=2)
  } else {
    title(main=main, cex.main=1)
  }

  ## Add ticks
  add_date_ticks(res)
  return(NULL)
}


##' For a matrix of CV scores (which are included in output from the function
##' \code{aggregateres()} or \code{blockcv_summary()}, make a 2d heatmap
##' @param cvscore.mat Matrix containing CV scores.
plot_cvscore <- function(cvscore.mat){
  mat = cvscore.mat
  colnames(mat) = signif(as.numeric(colnames(mat)),2)
  rownames(mat) = signif(as.numeric(rownames(mat)),2)
  drawmat_precise(mat, contour = FALSE,
                  ylab = expression(lambda[alpha]),
                  xlab = expression(lambda[beta]))
}


##' Another helper function to /precisely/ draw the entries of a matrix.
##' @param mat Matrix of interest.
##' @param contour If \code{TRUE}, draw a contour using
##'   \code{lattice::levelplot()}.
##' @param ... Other arguments to \code{lattice::levelplot()}.
##'
##' @return lattice object.
drawmat_precise <- function(mat, contour = FALSE, ...){

## Dummy data
## data <- matrix(runif(100, 0, 5) , 10 , 10)

  if(is.null(colnames(mat))){
    colnames(mat) <- paste(rep("col\n",ncol(mat)),
                            c(1:ncol(mat)) , sep=" ")
    rownames(mat) <- paste(rep("row",nrow(mat)),
                           c(1:nrow(mat)) , sep=" ")
  }

  ## Color function
  colfun = colorRampPalette(c("blue", "red"))

  # plot it flipping the axis
  lattice::levelplot(t(mat[c(nrow(mat):1) , ]),
                     col.regions = colfun(100),
                     contour = contour,
                     ## xaxt = 'n',
                     las = 2,
                     ...)
}


##' [MORE GENERIC VERSION] Only in the 2d case, make series of plots.
##'
##' @param show.fewer NULL by default. Otherwise the indices of the time points
##'   whose plots are to be shown.
##' @param obj Result of running covarem().
##' @param ylist List of responses.
##'
##' @return Nothing
plot2d_generic.covarem <- function(obj, ylist, ask=TRUE, show.fewer=NULL,
                                   show.xb.constraint=FALSE,
                                   dims=c(1:2)
                                   ## maxdev=NULL## Temporary addition.
                                   ){

  ## Basic checks
  ## assert_that(obj$dimdat==2)

  ## Collect plot ranges
  all.y = do.call(rbind, ylist)
  ylim = range(all.y[,2])
  xlim = range(all.y[,1])

  ## Extract numclust  (eventually from somewhere else)
  mns = obj$mn[,dims,]
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

  if(ask) par(ask=TRUE)
  ## par(mfrow=c(1,3))
  m = matrix(c(1, 1, 2, 2,
               1, 1, 3, 3),
             nrow = 2, ncol = 4, byrow=TRUE)
  ## m = matrix(c(1, 1, 2, 2,
  ##              1, 1, 2, 2,
  ##              3, 3, 4, 4,
  ##              3, 3, 5, 5),
  ##            nrow = 2, ncol = 4, byrow=TRUE)
  layout(m)
  par(mar=c(3,2,4,1))
  for(tt in TTrange){

    ## 1. Add main plot of data and fitted means

    ## Create empty plot
    main0 = paste0("time ", tt, " out of ", TT)
    plot(NA, ylim=ylim, xlim=xlim, cex=3, ylab="", xlab="", main=main0,  font.main = 1)

    ## Add datapoints
    points(ylist[[tt]][,dims], col='grey90', pch=16, cex=.5)

    ## Add the ball constraint
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
      lines(ellipse::ellipse(x = obj$sigma[ kk,,],
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
  }
  par(ask=FALSE)
}




##' Only in the 2d case, make series of plots.
##'
##' @param show.fewer NULL by default. Otherwise the indices of the time points
##'   whose plots are to be shown.
##' @param obj Result of running covarem().
##' @param ylist List of responses.
##'
##' @return No return.
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
      lines(ellipse::ellipse(x = obj$sigma[kk,,],
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
    legend("topright", col=1:p, lwd=lwd, lty=1, legend=c("SST", "Salinity",
                                                         "Iron",
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
                     "SST", "",
                      "Salinity", "",
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
                     "SST",
                      "Salinity",
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
                     "SST", "",
                      "Salinity", "",
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
                     "SST",
                      "Salinity",
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






## ##' Default printing functionality for covarem() results.
## ##'
## ##' @param obj An object of class "covarem"
## ##'
## ##' @return No return.
## plot.covarem <- function(obj){

##   ## Setup
##   betas = do.call(cbind, lapply(obj$beta, function(beta)beta[-1,]))
##   alpha =obj$alpha[ ,-1]
##   cols = RColorBrewer::brewer.pal(obj$numclust, "Set3")
##   dimdat = obj$dimdat

##   ## Plot means, in each dimension
##   par(ask=TRUE)
##   for(idim in 1:dimdat){
##     mn = obj$mn[,idim,]
##     par(mfrow=c(1,2))
##     matplot(mn, type='l', lty=1, lwd=2, main=paste0("Means in ", idim, "-th dim"),
##             col=cols)
##     wt = obj$pie
##     matplot(wt, type='l', lty=1, lwd=2, main="Weight", ylim=c(0,1),
##             col=cols)
##   }

##   ## Plot objective values
##   plot(res$obj[-1], type='o', main="Objectives", pch=16, ylab="value", xlab="iterations")
##   par(ask=FALSE)
## }
