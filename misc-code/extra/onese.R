
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
onese <- function(destin,
                  gridsize = 12){
## outputdir="~/repos/flowcy/output/vanilla"

  ## Gather the degrees of freedom
  dfmat = aggregateres_df(gridsize, destin)
  cvscoremat = aggregateres(gridsize, destin)

  ## Add the line
  ind = ind.min = which(cvscoremat == min(cvscoremat, na.rm=TRUE), arr.ind=TRUE)

  ## Add horizontal lines
  load(file=file.path(destin, paste0(ind[1], "-", ind[2], ".Rdata")))
  se = sd(cvres$all)/sqrt(length(cvres$all))

  ## Apply the 1SE rule
  inds = which((cvscoremat <  cvscoremat[ind[1], ind[2]] + se &
                cvscoremat >  cvscoremat[ind[1], ind[2]] - se &
               dfmat < dfmat[ind[1], ind[2]]),
               arr.ind = TRUE)
  if(nrow(inds) > 0){

    ## Pick the minimum CV error that is most parsimonious.
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
              res = res,
              se = se,
              ind.min = ind.min,
              cvscoremat = cvscoremat,
              dfmat = dfmat))
}
