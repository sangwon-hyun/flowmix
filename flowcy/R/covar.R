##' Main function for covariate EM.
##' @param ylist  T-length list each containing response matrices  of size (nt x
##'   3), which contains coordinates of  the 3-variate particles, organized over
##'   time (T) and with (nt) particles at every time.
##' @param X Matrix of size (T x p+1)
##' @param pie.list (T by K)
##' @param mean_lambda lambda for lasso for the mean.
##' @param pie_lambda lambda for lasso for pie.
##' @return  List containing  fitted parameters and  means and  mixture weights,
##'   across algorithm iterations.
covarem <- function(ylist, X=NULL, numclust, niter=100, mn=NULL, pie_lambda=0,
                    mean_lambda=0, verbose=FALSE,
                    warmstart = c("none", "rough"), sigma.fac=1, tol=1E-6){

  ## Setup.
  ntlist = sapply(ylist, nrow)
  dimdat = ncol(ylist[[1]])
  TT = length(ylist)
  p = ncol(X)
  warmstart = match.arg(warmstart)

  ## Initialize.
  beta = init_beta(TT, p, dimdat, numclust)
  alpha = init_alpha(dimdat, p)
  if(is.null(mn)){
    if(warmstart=="rough"){
      mn = warmstart_covar(ylist, numclust)
    } else if (warmstart=="none"){
      mn = aperm(init_mu(lapply(ylist, cbind), numclust, TT), c(1,3,2))
    } else {
      stop("warmstart option not recognized")
    }
  }
  pie = calc_pie(TT, numclust) ## Let's just say it is all 1/K for now.
  sigma = init_sigma(ylist, numclust, TT, fac=sigma.fac) ## (T x numclust x dimdat x dimdat)

  ## Initialize alpha and beta
  beta.list = alpha.list = sigma.list = pie.list = mn.list = list()
  objectives = rep(NA, niter)
  objectives[1] = -1E20 ## Fake
  beta.list[[1]] = beta ## beta.list: Each element is a (T x p+1 x 3 x K) array
  alpha.list[[1]] = alpha ## alpha.list: Each element is a T by p+1 array
  mn.list[[1]] = mn
  sigma.list[[1]] = sigma
  pie.list[[1]] = pie

  start.time=Sys.time()
  for(iter in 2:niter){
    if(verbose) printprogress(iter, niter, "EM iterations.", start.time=start.time)

    ## Conduct E step
    resp <- Estep_covar(mn.list[[iter-1]],
                        sigma.list[[iter-1]],
                        pie.list[[iter-1]],
                        ylist,
                        numclust,
                        ntlist)  ## This should be (T x numclust x dimdat x dimdat)

    ## Conduct M step
    ## 1. Alpha
    res.alpha = Mstep_alpha(resp,
                            X, numclust,
                            lambda=pie_lambda)
    alpha.list[[iter]] = res.alpha$alpha
    pie.list[[iter]] = res.alpha$pie

    ## 2. Beta
    res.beta = Mstep_beta_faster_lasso(resp, ylist, X,
                                       mean_lambda=mean_lambda,
                                       sigma.list[[iter-1]])
    beta.list[[iter]] = res.beta$beta
    mn.list[[iter]]    = res.beta$mns

    ## 3. Sigma
    sigma.list[[iter]] <- Mstep_sigma_covar(resp,
                                            ylist,
                                            mn.list[[iter]],
                                            numclust)

    ## Calculate the objectives
    objectives[iter] = objective_overall_cov(aperm(mn.list[[iter]], c(1,3,2)),
                                             pie.list[[iter]],
                                             sigma.list[[iter]],
                                             ylist,
                                             pie_lambda=pie_lambda,
                                             mean_lambda=mean_lambda,
                                             alpha=res.alpha$alpha,
                                             beta=res.beta$beta)

    ## Check convergence
    if(check_converge_rel(objectives[iter-1],
                          objectives[iter], tol=tol)) break
  }

  ## Threshold (now that we're using CVXR for beta)
  beta=beta.list[[iter]]
  betathresh = 1E-3
  beta = lapply(beta, function(a){
    a[abs(a)<betathresh] = 0
    a
  })

  return(list(alpha=alpha.list[[iter]],
              alpha.fit=res.alpha$fit,
              beta=beta,
              mn=mn.list[[iter]],
              pie=pie.list[[iter]],
              sigma=sigma.list[[iter]],
              objectives=objectives[1:iter],
              final.iter=iter,
              ## Above is output, below are data/algorithm settings.
              ntlist=ntlist,
              dimdat=dimdat,
              TT=TT,
              p=p,
              numclust=numclust,
              X=X,
              pie_lambda=pie_lambda,
              mean_lambda=mean_lambda
              ))
}


##' Prediction: given  new X's,  generate a set of means and pies (and return
##' the same Sigma)
##' @param res object returned from covariate EM covarem().
predict.covarem <- function(res, newx=NULL){

  ## ## Check the dimensions
  ## stopifnot(ncol(new.x) == ncol(res$X))
  ## newx = X[1,,drop=FALSE]
  if(is.null(newx)){
    newx = res$X
  }

  ## Augment it with a dummy variable 1
  if(nrow(newx)>1){
    newx.a = cbind(rep(1, nrow(newx)), newx)
  } else {
    newx.a = c(1, newx)
  }

  ## Predict the means (manually).
  newmn = lapply(1:numclust, function(iclust){
    newx.a  %*%  res$beta[[iclust]]
  })
  newmn = abind::abind(newmn, along=0)
  newmn = aperm(newmn, c(2,3,1)) ## This needs to by (T x dimdat x numclust)

  ## Predict the pies.
  newpie = predict(res$alpha.fit, newx=newx, type='response')[,,1]

  ## Return all three things
  return(list(newmn=newmn,
              newpie=newpie,
              sigma=res$sigma))
}

##' Plot results from covariate EM \code{covarem()}.
plot.covarem <- function(res, newx=NULL){

  ## Extract numclust  (eventually from somewhere else)
  mns = res$mn.list[[plot.iter]]
  numclust = res$numclust

  ## Define range of date
  xlim = range(do.call(rbind, res$ylist)[,1])
  ylim = range(do.call(rbind, res$ylist)[,2])

  ## General plot settings
  cols = RColorBrewer::brewer.pal(numclust, "Set3")

  ## Make /series/ of plots.
  par(ask=TRUE)
  par(mfrow=c(1,3))
  for(tt in 1:res$TT){
    main0 = paste0("time ", tt, " out of ", TT)
    main = paste0("Iteration ", plot.iter)
    plot(NA, ylim=ylim, xlim=xlim, cex=3, ylab="", xlab="", main=main0)

    ## Add datapoints
    points(res$ylist[[tt]], col='lightgrey', pch=16, cex=.5)

    ## Add truth as well
    if(!is.null(truths)){
      numclust.truth = dim(truths$mns)[3]
      for(kk in 1:numclust.truth){
        points(x=truths$mns[tt,1,kk],
               y=truths$mns[tt,2,kk],
               col="black", pch=16, cex=truths$pies[[kk]][tt]*5)
      }
    }
  }

}

## ##' Visualizes the  data and  results of  running mixture  EM. REQUIRES  a |res|
## ##' object.
## make_2d_plot <- function(figdir, plotname, res, plot.iter=NULL, truths=NULL, h=8){

##   if(is.null(plot.iter)) plot.iter = res$final.iter
##   ## pdf(file=file.path(figdir, plotname), width=15, height=8)
##   pdf(file=file.path(figdir, plotname), width=3*h, height=h)

##   ## Extract numclust  (eventually from somewhere else)
##   mns = res$mn.list[[plot.iter]]
##   numclust = res$numclust


##   ## Define range of date
##   xlim = range(do.call(rbind, res$ylist)[,1])
##   ylim = range(do.call(rbind, res$ylist)[,2])

##   ## General plot settings
##   cols = RColorBrewer::brewer.pal(numclust, "Set3")

##   ## Eventually extract TT some other way.
##   par(mfrow=c(1,3))
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

##     ## Collect pies
##     pies = lapply(1:numclust, function(iclust){
##       res$pie.list[[plot.iter]][,iclust]
##     })

##     ## Add fitted means
##     for(kk in 1:numclust){
##       points(x=mns[tt,1,kk],
##              y=mns[tt,2,kk],
##              col='red', pch=16, cex=pies[[kk]][tt]*5)
##     }

##     ## Add legend
##     legend("bottomright", col=c("black", "red"),
##            pch=c(16,16), legend = c("truth", "fitted"))

##     for(kk in 1:numclust){
##     lines(ellipse::ellipse(x=res$sigma.list[[plot.iter]][tt,kk,,],
##                            centre=mns[tt,,kk]
##                            ), lwd=1/2, col='red')
##     }
##     ## lines(ellipse(rho), col="red")       # ellipse() from ellipse package
##     ## lines(ellipse(rho, level = .99), col="green")
##     ## lines(ellipse(rho, level = .90), col="blue")


##     ## Plot pies
##     plot(NA, xlim=c(0,TT), ylim=c(0,1.4), main=main,
##          ylab = "Mixture probability", xlab="time, t=1,..,T")
##     for(iclust in 1:numclust){
##       lines(pies[[iclust]], col=cols[iclust], lwd=2, type='l')
##     }
##     abline(v=tt, col='green')


##     plot(NA, xlim=c(0,TT), ylim=range(res$X)*1.5, main=main,
##          ylab = "Environmental covariates", xlab="time, t=1,..,T")
##     for(ii in 1:2){
##       lines(res$X[,ii], col=ii, lwd=2, type='l')
##     }
##     abline(v=tt, col='green')

##     ## if(plot.iter>1){
##     ##   betas = res$beta.list[[plot.iter]]
##     ##   betas = round(do.call(cbind, betas),2)
##     ##   rownames(betas)[1] = "intp"
##     ##   colnames(betas)= paste0("cluster ", c(1,1,2,2))
##     ##   text(0.6,1, paste(capture.output(betas), collapse='\n'), pos=4)##, family="monospace")

##     ##   alphas = as.matrix(t(res$alpha.list[[plot.iter]]))
##     ##   rownames(alphas)[1] = "intp"
##     ##   alphas = round(alphas,2)
##     ##   colnames(alphas)= paste0("cluster ", c(1,2))
##     ##   text(35,1, paste(capture.output(alphas), collapse='\n'), pos=4)##, family="monospace")
##     ## }
##   }
##   graphics.off()
## }



## ##' Visualizes the  data only
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
