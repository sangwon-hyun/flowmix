## Synopsis: contains the main CV wrapper, and all helper functions related to
## cross-validations.

##' CV wrapper for covarem().
##' @param nsplit Number of CV splits. Defaults to 5.
##' @param ... default arguments to covarem().
##' @return List containing (1) the set of coefficients
cv.covarem <- function(ylist=ylist, X=X, mean_lambdas=NULL,
                       pie_lambdas=NULL,
                       max_mean_lambda=NULL,
                       max_pie_lambda=NULL,
                       gridsize=5,
                       nsplit=5,
                       verbose=FALSE,
                       refit=FALSE,
                       multicore.cv=FALSE,
                       ...){
  ## Basic checks
  stopifnot(length(mean_lambdas) == length(pie_lambdas))
  if(is.null(mean_lambdas)){
    assert_that(!is.null(max_mean_lambda))
    mean_lambdas = c(exp(seq(from=-8, to=log(max_mean_lambda), length=gridsize)))
  }
  if(is.null(pie_lambdas)){
    assert_that(!is.null(max_pie_lambda))
    pie_lambdas = c(exp(seq(from=-8, to=log(max_pie_lambda), length=gridsize)))
  }

  ## Create CV split indices
  assert_that(nsplit >= 2)
  mysplits = cvsplit(ylist, nsplit=nsplit)

  ## Empty containers
  meanmat = sdmat = matrix(NA, ncol=gridsize, nrow=gridsize)
  cvscoremat = array(NA, dim=c(gridsize,gridsize, nsplit))

  reslist = list()
  start.time=Sys.time()
  for(ii in gridsize:1){
    for(jj in gridsize:1){
      if(verbose){
        print(c(ii,jj))
        print(Sys.time() - start.time)
        print("progressed")
      }

      ## Set warm starts.
      if(jj == gridsize){
        if(ii == gridsize)  warmstart_mn = NULL
        if(ii != gridsize)  warmstart_mn = reslist[[paste0(ii+1, "-", jj)]]$mn
      } else {
        warmstart_mn = reslist[[paste0(ii, "-", jj+1)]]$mn
      }

      ## Do CV to get score
      cvres = get_cv_score(ylist, X, mysplits, nsplit,
                           mean_lambda=mean_lambdas[ii],
                           pie_lambda=pie_lambdas[jj],
                           mn=warmstart_mn,
                           refit=refit,
                           multicore.cv=multicore.cv,
                           ...)

      ## Get the fitted results on the entire data
      res = covarem(ylist=ylist, X=X,
                    mean_lambda=mean_lambdas[ii],
                    pie_lambda=pie_lambdas[jj],
                    mn=warmstart_mn, ...)

      ## Store the results
      reslist[[paste0(ii, "-", jj)]] = res

      ## Record CV score
      meanmat[ii, jj] = cvres$mean
      sdmat[ii, jj] = sd(cvres$all)
      cvscoremat[ii,jj,]=cvres$all
    }
  }
  return(list(reslist=reslist,
              meanmat=meanmat,
              sdmat=sdmat,
              cvscoremat=cvscoremat))
}

##' Makes all cross-validation cytogram indices, by splitting the code{TT}
##' cytograms. Five groups of datasets are made according to the cytograms
##' holding the following indices: $\{(1,6,11,..), (2,7,12,...)  ...,
##' (5,10,15,...)\}$. Then for cross validation, pick one "test" group out of
##' the five, fit the data on the 4 other groups combined ("training" group),
##' and assess the likelihood of the model on the test group.
##' @param ylist ylist.
##' @param nsplit Defaults to 5.
##' @return \code{nsplit}-lengthed list, each containing groups of indices for
##'   cytograms.
cvsplit2 <- function(ylist, nsplit=5){
  TT = 100
  lapply(1:nsplit, function(isplit){
    (1:(TT/nsplit)) * nsplit + isplit
  })
}


##' Makes all cross-validation indices, by splitting each cytogram's \eqn{n_t}
##' particles 5 ways, and storing the 5 splitted indices into a TT-lengthed
##' list.
##' @param ylist ylist.
##' @param nsplit Defaults to 5.
##' @return \code{TT}-lengthed list, each containing 5 groups of indices split
##'   \code{nsplit}-ways -- these are the folds to be used during cross
##'   validation.
cvsplit<- function(ylist, nsplit = 5){
  TT = length(ylist)
  all.cv.inds = lapply(1:TT, function(tt){
    flds <- caret::createFolds(ylist[[tt]][,1], k = nsplit, list = TRUE, returnTrain = FALSE)
    names(flds)[1] <- "train"
    return(flds)
  })
}

##' Getting cross-validated test likelihood, averaged over train/test splits.
##' @param splits TT-lengthed list of indices.
##' @param ylist List of responses.
##' @param ... arguments to covarem.
##' @return Cross validated test likelihood, scalar-valued.
get_cv_score <- function(ylist, X, splits, nsplit, refit,
                         multicore.cv = FALSE,
                         ...){
  ## stopifnot(length(splits[[1]])!=nsplit) ## good check but only works if TT>1
  if(multicore.cv){
    mc.cores = nsplit
  } else {
    mc.cores = 1
  }

  ## Cycle through splits, and calculate CV scroe
  all.scores = mclapply(1:nsplit, function(test.isplit){
    get_cv_score_onesplit(test.isplit, splits, ylist, X, refit,## , refit, sel_coef,
                          ...)
  }, mc.cores=mc.cores)
  all.scores = do.call(c, all.scores)
  return(list(mean=mean(all.scores), all=all.scores))
}

##' Inner  function  for calculating  cross-validated  test  likelihood for  one
##' train/test split.   Specifically, it takes  the set of splitted  indices for
##' all times in 1:TT, trains on (1:nsplit)[-isplit], tests on splits
##' @param splits TT-lengthed list of indices.
##' @param test.split The split (out of 1:nsplit) to use for test.
##' @param ylist List of responses.
##' @param refit (experimental), defaults to FALSE. If TRUE, then the refitted
##'   non-regularized solutions (with only a user-specified set of active
##'   coefficients) are calculated.
##' @param ... arguments to covarem
##' @return One split's test likelihood.
get_cv_score_onesplit <- function(test.isplit, splits, ylist, X, refit, ...){##, refit=FALSE,...){

  TT = length(ylist)
  dimdat = ncol(ylist[[1]])
  myargs = list(...)
  numclust = myargs$numclust ## might work, might not.

  ## Obtain train and test data according to test split
  ylist.test = lapply(1:TT, function(tt){
    ind = splits[[tt]][[test.isplit]]
    ylist[[tt]][ind,]
  })
  ylist.train = lapply(1:TT, function(tt){
    ind = splits[[tt]][[test.isplit]]
    ylist[[tt]][-ind,]
  })

  ## Run algorithm on training data, evaluate on test data.
  res.train = covarem(ylist.train, X,  refit = FALSE, ...)

  ## If applicable, unregularized refit on the sparsity pattern
  if(refit){
    sel_coef = get_sparsity_pattern(res.train)
    res.train = covarem(ylist.train, X,  refit = TRUE,
                        sel_coef = sel_coef, ...)
  } else {
    sel_coef = NULL
  }

  ## Assign mn and pie
  pred = predict.covarem(res.train)
  stopifnot(all(pred$newpie >= 0))

  ## Calculate objective (penalized likelihood)
  objective_overall_cov(pred$newmn,
                        pred$newpie,
                        pred$sigma,
                        TT, dimdat, numclust,
                        ylist.test,
                        pie_lambda = 0,
                        mean_lambda = 0,
                        alpha = res.train$alpha,
                        beta = res.train$beta)
}


##' Create a list containing the candidate regularization parameters for EM.
make_lambdas <- function(ylist, X, numclust, cv.grid.size=5){
  res0 = get_param0(ylist, numclust)
  max_lambda_beta = lambda_beta_max(res0$alpha0, res0$beta0, res0$sigmalist0,
                                    numclust, X, res0$resplist0, ylist)
  lambdas_beta = c(0,exp(seq(from=0, to=log(max_lambda_beta),
                             length=cv.grid.size)))
  max_lambda_alpha = lambda_alpha_max(res0$alpha0, res0$beta0, res0$sigmalist0,
                                      numclust, X, ylist)
  lambdas_alpha = c(0,exp(seq(from=0, to=log(max_lambda_alpha), length=cv.grid.size)))
  return(list(lambdas_beta=lambdas_beta,
              lambdas_alpha=lambdas_alpha))
}


##' Helper to get sparsity pattern of fitted coefficients.
##' @param res covarem class object.
##' @return list containing two objects (alpha and beta) that are the same
##'   structure as \code{res$alpha} and \code{beta}, but are boolean
##'   matrices. Entries are TRUE if they are to be selected, and FALSE they are
##'   zero.
get_sparsity_pattern <- function(res, thresh=1E-8){

  list(beta = lapply(res$beta, function(onebeta){as.matrix(abs(onebeta) <= thresh)}),
       alpha = as.matrix(abs(res$alpha) <= thresh))
}



##' Make a plot of the 1SE results.
onese_plot <- function(isim, outputdir,
                       gridsize=12){

  obj = onese(isim, outputdir, gridsize)
  ind.min = obj$ind.min
  dfmat = obj$dfmat
  cvscoremat = obj$cvscoremat
  se = obj$se

  ## Initialize the plot
  plot(NA,
       xlim=c(0, max(dfmat)),
       ylim=c(min(cvscoremat, na.rm=TRUE), max(cvscoremat, na.rm=TRUE))+c(-1000,1000),
       ylab="CV error",
       xlab="Total # nonzero coefs")

  ## Add the points
  xlist = c()
  ylist = c()
  for(ii in 1:gridsize){
    for(jj in 1:gridsize){
      xlist = c(xlist, dfmat[ii,jj])
      ylist = c(ylist, cvscoremat[ii,jj])
      points(x=dfmat[ii,jj],
             y=cvscoremat[ii,jj])
    }
  }

  ## if(any(is.na(cvscoremat))) cvscoremat[is.na(cvscoremat)]=Inf
  ind = which(cvscoremat==min(cvscoremat, na.rm=TRUE), arr.ind=TRUE)
  abline(v=dfmat[ind[1], ind[2]], col='green')
  abline(h=cvscoremat[ind[1], ind[2]] + c(-se, se), col='violet')

  points(x = dfmat[ind.min[1], ind.min[2]],
         y = cvscoremat[ind.min[1], ind.min[2]], pch = 16,
         col = 'violet', cex = 1.5)

  ## Add kernel density plot
  if(any(is.na(ylist))){ xlist = xlist[-which(is.na(ylist))]; ylist = ylist[-which(is.na(ylist))]}
  a = ksmooth(x=xlist, y=ylist, kernel = "normal", bandwidth=5)
  lines(a$y~a$x, col='red')
}

##' Apply the 1SE rule (results are similar to those of get_optimal_info())
##' @param outputdir Location of the output files. e.g. outputdir="~/output/vanilla"
##' @return NULL.
onese <- function(isim, outputdir,
                       gridsize=12){
## outputdir="~/repos/flowcy/output/vanilla"

  ## Load data
  destin = file.path(outputdir, paste0("isim-", isim))

  ## Gather the degrees of freedom
  dfmat = aggregateres_df(gridsize, destin)
  cvscoremat = aggregateres(gridsize, destin)

  ## Add the line
  ind = ind.min = which(cvscoremat==min(cvscoremat, na.rm=TRUE), arr.ind=TRUE)

  ## Add horizontal lines
  load(file=file.path(destin, paste0(ind[1], "-", ind[2], ".Rdata")))
  se = sd(cvres$all)/sqrt(length(cvres$all))

  ## Apply the 1SE rule
  inds = which((cvscoremat <  cvscoremat[ind[1], ind[2]] + se &
                cvscoremat >  cvscoremat[ind[1], ind[2]] - se &
               dfmat < dfmat[ind[1], ind[2]]),
               arr.ind = TRUE)
  if(nrow(inds) > 0){

    ## Pick the minimum CV error that is most sparsimonious.
    all.df = apply(inds, 1, function(myind){dfmat[myind[1],myind[2]]})
    inds = inds[which(all.df==min(all.df)),,drop=FALSE]
    all.cvscores = apply(inds, 1, function(myind){cvscoremat[myind[1],myind[2]]})
    ind.min = inds[which.min(all.cvscores),,drop=FALSE]
    assert_that(nrow(ind.min)==1)

  }

  ## Obtain the best parameters and return
  filename = paste0(ind.min[1], "-", ind.min[2], ".Rdata")
  load(file=file.path(destin, filename))
  lambda_alpha = alpha_lambdas[ind.min[1]]
  lambda_beta = beta_lambdas[ind.min[2]]

  return(list(param = c(lambda_alpha, lambda_beta),
              res=res,
              se=se,
              ind.min=ind.min,
              cvscoremat=cvscoremat,
              dfmat=dfmat))

}



##' Makes a 'fancy' plotly plot.
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
