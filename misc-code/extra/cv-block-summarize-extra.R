##' Apply the 1SE rule (results are similar to those of get_optimal_info())
##' @param outputdir Location of the output files. e.g. outputdir="~/output/vanilla"
##' @return NULL.
blockcv_onese <- function(destin,
                          gridsize,
                          nrep,
                          sim=FALSE, isim=NULL){


  ## Pseudocode
  ## 1. Get CV score matrix.
  ## 2. Get the DFs.
  ## 3. Get the SD of all the folds AT the bestres.
  ## 4. Take all the models that are within 1 SE of the bestres.
  ## 5. Calculate the DFs of these models.
  ## 6. Find the most parsimonious model this way.

  ## 1. Get CV score matrix.
  destin = "~/Dropbox/research/usc/hpc-output/blockcv-2-75-5"
  nfold = 5
  nrep = 5
  cv_gridsize=10
  obj = blockcv_aggregate(destin=destin, cv_gridsize=cv_gridsize, nfold=nfold, nrep=nrep, sim = FALSE)
  ialpha = obj$min.inds[1]
  ibeta = obj$min.inds[2]
  cvscore.mat = obj$cvscore.mat

  ## 2. Get the DFs
  dfmat = blockcv_aggregate_df(gridsize=cv_gridsize, nrep=nrep, destin=destin)$mat

  ##3. Get the SD of all the folds at the min.ind
  dim(obj$cvscore.array)
  se = sd(apply(obj$cvscore.array[min.inds[1], min.inds[2],,],1,min))

   ## 4. Get all the models within 1SE of the bestres
  inds = which((cvscore.mat <  cvscore.mat[min.inds[1], min.inds[2]] + se &
                cvscore.mat >  cvscore.mat[min.inds[1], min.inds[2]] - se &
                dfmat < dfmat[min.inds[1], min.inds[2]]),
               arr.ind = TRUE)

  ## 5. Get all the DFs of these models
  if(nrow(inds) > 0){

    ## Pick the minimum CV error that is most parsimonious.
    all.df = apply(inds, 1, function(myind){dfmat[myind[1],myind[2]]})
    inds = inds[which(all.df==min(all.df)),,drop=FALSE]
    all.cvscores = apply(inds, 1, function(myind){cvscore.mat[myind[1],myind[2]]})
    onese.min.inds = inds[which.min(all.cvscores),,drop=FALSE] ## Out of the lowest DF ties, get the smallest CV score.
    assert_that(nrow(onese.min.inds)==1)

  }

  ## ## Gather the degrees of freedom
  ## ## dfmat = aggregateres_df(gridsize, destin)
  ## ## cvscoremat = aggregateres(gridsize, destin)
  ## cvscore.mat = blockcv_aggregate(destin = destin,
  ##                                 cv_gridsize = cv_gridsize,
  ##                                 nfold = nfold,
  ##                                 sim = TRUE,
  ##                                 isim = 1,
  ##                                 nrep = nrep)$cvscore.mat

  ## ## Add the line
  ## ind = ind.min = which(cvscoremat == min(cvscoremat, na.rm=TRUE), arr.ind=TRUE)

  ## ## Add horizontal lines
  ## load(file=file.path(destin, paste0(ind[1], "-", ind[2], ".Rdata")))
  ## se = sd(cvres$all)/sqrt(length(cvres$all))

  ## ## Apply the 1SE rule
  ## inds = which((cvscoremat <  cvscoremat[ind[1], ind[2]] + se &
  ##               cvscoremat >  cvscoremat[ind[1], ind[2]] - se &
  ##               dfmat < dfmat[ind[1], ind[2]]),
  ##              arr.ind = TRUE)
  ## if(nrow(inds) > 0){

  ##   ## Pick the minimum CV error that is most parsimonious.
  ##   all.df = apply(inds, 1, function(myind){dfmat[myind[1],myind[2]]})
  ##   inds = inds[which(all.df==min(all.df)),,drop=FALSE]
  ##   all.cvscores = apply(inds, 1, function(myind){cvscoremat[myind[1],myind[2]]})
  ##   ind.min = inds[which.min(all.cvscores),,drop=FALSE]
  ##   assert_that(nrow(ind.min)==1)

  ## }

  ## Obtain the best parameters and return
  ## filename = paste0(ind.min[1], "-", ind.min[2], ".Rdata")
  ## load(file=file.path(destin, filename))
  lambda_alpha = lambda_alphas[onese.min.inds[1]]
  lambda_beta = lambda_alphas[onese.min.inds[2]]

  ## Obtain the best model
  load(file.path(destin, "bestreslist.Rdata"))## If necessary, load or calculate the bestreslist
  bestres_onese = bestreslist[[paste0(onese.min.inds[1], "-", onese.min.inds[2])]]

  return(list(param = c(lambda_alpha, lambda_beta),
              bestres_onese = bestres_onese,
              se = se,
              min.inds = min.inds,
              onese.min.inds = onese.min.inds,
              inds = inds,
              cvscoremat = cvscoremat,
              dfmat = dfmat))
}
