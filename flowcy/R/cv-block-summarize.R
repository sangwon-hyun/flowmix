##' Get rid of the |nrep| replicates by retaining only the best one.
##' 1. If there are only nrep guys we will
blockcv_reduce_by_nrep <- function(destin, cv_gridsize, nfold, nrep,
                              sim=FALSE, isim = NULL,
                              save=FALSE, resfile = "all-cvres.Rdata"){

  ## Read the meta data (for |nfold|, |cv_gridsize|, |nrep|)
  load(file = file.path(destin, 'meta.Rdata'))

  ## Identify the best replicate
  cvscore.array = array(NA, dim = c(cv_gridsize, cv_gridsize, nfold, nrep))
  cvscore.mat = matrix(NA, nrow = cv_gridsize, ncol = cv_gridsize)
  for(ialpha in 1:cv_gridsize){
    for(ibeta in 1:cv_gridsize){
      ## cvscores <- sapply(1:nfold, function(igrid){
      obj = matrix(NA, nrow=nfold, ncol=nrep)
      for(ifold in 1:nfold){
        for(irep in 1:nrep){
          filename = make_cvscore_filename(ialpha, ibeta, ifold, irep, sim, isim)
          tryCatch({
            load(file.path(destin, filename))
            cvscore.array[ialpha, ibeta, ifold, irep] = cvscore
            obj[ifold, irep] = objectives[length(objectives)]
          }, error = function(e){})

          ## ## Also compare to the aggregated result so far.
          ## aggregate_filename = make_cvscore_filename(ialpha, ibeta, ifold, irep = NULL, sim, isim)
          ## ## if irep=NULL, look for "aggr-" prefix; bake this into the make_cvscore_filename() as well.
          ## load(file.path(destin, filename))
          ## obj_aggr = objectives[length(objectives)]
        }
      }

      ## Pick out the CV scores with the *best* (lowest) objective value
      cvscores = cvscore.array[ialpha, ibeta,,]
      best.models = apply(obj, 1, function(myrow) which(myrow == min(myrow, na.rm=TRUE)))
      final.cvscores = sapply(1:nfold, function(ifold){
        cvscores[ifold, best.models[ifold]]
      })
      cvscore.mat[ialpha, ibeta] = mean(final.cvscores)
    }
  }

  ## Clean a bit
  cvscore.mat[which(is.nan(cvscore.mat), arr.ind=TRUE)] = NA
  rownames(cvscore.mat) = signif(pie_lambdas,3)
  colnames(cvscore.mat) = signif(mean_lambdas,3)


  ## Find the minimum
  mat = cvscore.mat
  min.inds = which(mat == min(mat, na.rm=TRUE), arr.ind=TRUE)

  ## Recent addition
  if(save){
    cat("Saving aggregated results to ", file.path(destin, resfile))
    save(cvscore.array,
         cvscore.mat,
         mean_lambdas,
         pie_lambdas,
         min.inds,
         file = file.path(destin, resfile))
  }

  return(list(cvscore.array = cvscore.array,
              cvscore.mat = cvscore.mat,
              mean_lambdas = mean_lambdas,
              pie_lambdas = pie_lambdas,
              min.inds = min.inds))
}



##' Aggregate results from blocked CV.
##'
##' @param destin destination folder for output.
##' @param cv_gridsize Grid size for cross validation.
##' @param nfold Number of folds.
##'
##' @export
## blockcv_aggregate <- function(destin, cv_gridsize, nfold){
blockcv_aggregate <- function(destin, cv_gridsize, nfold, nrep,
                              sim=FALSE, isim = NULL,
                              save=FALSE, resfile = "all-cvres.Rdata"){

  ## Read the meta data (for |nfold|, |cv_gridsize|, |nrep|)
  load(file = file.path(destin, 'meta.Rdata'))

  ## Aggregate the results
  cvscore.array = array(NA, dim = c(cv_gridsize, cv_gridsize, nfold, nrep))
  cvscore.mat = matrix(NA, nrow = cv_gridsize, ncol = cv_gridsize)
  for(ialpha in 1:cv_gridsize){
    for(ibeta in 1:cv_gridsize){
      ## cvscores <- sapply(1:nfold, function(igrid){
      obj = matrix(NA, nrow=nfold, ncol=nrep)
      for(ifold in 1:nfold){
        for(irep in 1:nrep){
          ## if(ialpha == 2 & ibeta == 1 & ifold ==5 & irep == 5)
          ## ialpha = 2 ; ibeta = 1 ; ifold =5 ; irep = 5
          ## filename = paste0(ialpha, "-", ibeta, "-", ifold,
          ##                   "-", irep, "-cvscore.Rdata")
          filename = make_cvscore_filename(ialpha, ibeta, ifold, irep, sim, isim)
          tryCatch({
            load(file.path(destin, filename))
            ## print(filename)
            cvscore.array[ialpha, ibeta, ifold, irep] = cvscore
            obj[ifold, irep] = objectives[length(objectives)]
          }, error = function(e){})
        }
      }

      ## Pick out the CV scores with the *best* (lowest) objective value
      cvscores = cvscore.array[ialpha, ibeta,,]
      best.models = apply(obj, 1, function(myrow){
        ind = which(myrow == min(myrow, na.rm=TRUE))
        if(length(ind)>1) ind = ind[1]  ## Just choose one, if there is a tie.
        return(ind)
      })
      final.cvscores = sapply(1:nfold, function(ifold){
        cvscores[ifold, best.models[ifold]]
      })
      cvscore.mat[ialpha, ibeta] = mean(final.cvscores)
    }
  }

  ## Clean a bit
  cvscore.mat[which(is.nan(cvscore.mat), arr.ind=TRUE)] = NA
  rownames(cvscore.mat) = signif(pie_lambdas,3)
  colnames(cvscore.mat) = signif(mean_lambdas,3)


  ## Find the minimum
  mat = cvscore.mat
  min.inds = which(mat == min(mat, na.rm=TRUE), arr.ind=TRUE)

  ## Recent addition
  if(save){
    cat("Saving aggregated results to ", file.path(destin, resfile), fill=TRUE)
    save(cvscore.array,
         cvscore.mat,
         mean_lambdas,
         pie_lambdas,
         min.inds,
         file = file.path(destin, resfile))
  }

  return(list(cvscore.array = cvscore.array,
              cvscore.mat = cvscore.mat,
              mean_lambdas = mean_lambdas,
              pie_lambdas = pie_lambdas,
              min.inds = min.inds))
}


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


##' Helper to aggregate parallelized CV results and obtain degrees of freedom
##' (DF) estimate, saved in |destin|.
##'
##' @param gridsize Size of CV grid.
##' @param destin Location of saved things.
##'
##' @return Matrix containing estimated degrees of freedom.
blockcv_aggregate_df <- function(gridsize, nrep, destin,
                                 save=FALSE, resfile = "all-cvres-df.Rdata"){

  df.array = obj.array = df.alpha.array = df.beta.array = array(NA, dim=c(gridsize, gridsize, nrep))
  df.mat = df.alpha.mat = df.beta.mat = matrix(NA, ncol=gridsize, nrow=gridsize)
  for(ialpha in 1:gridsize){
    for(ibeta in 1:gridsize){

      ## Objective value
      obj = rep(NA, nrep)
      df = df.alpha = df.beta = rep(NA, nrep)
      for(irep in 1:nrep){
        ## print(c(ialpha, "-", ibeta, "-", irep))

        tryCatch({
          ## Load fitted result
          filename = paste0(ialpha, "-", ibeta, "-", irep, "-fit.Rdata")
          load(file.path(destin, filename))

          ## Calculate DF
          df[irep] = do.call(sum, lapply(res$beta, function(mybeta){
            sum(mybeta[-1,]!=0)})) + sum(res$alpha[,-1]!=0)

          df.alpha[irep] = sum(res$alpha[,-1]!=0)

          df.beta[irep] = do.call(sum, lapply(res$beta, function(mybeta){
            sum(mybeta[-1,]!=0)}))

          ## Also calculate objective function
          objectives = res$objectives
          obj[irep] = objectives[length(objectives)]

        }, error = function(e){})
      }
      df.array[ialpha, ibeta, ] = df
      obj.array[ialpha, ibeta,] = obj
      df.alpha.array[ialpha, ibeta, ] = df.alpha
      df.beta.array[ialpha, ibeta, ] = df.beta

      ## Calculate the df of the best model
      if(!all(is.na(obj))){
        ## df.mat[ialpha, ibeta] = df[which.max(obj, na.rm=TRUE)]
        df.mat[ialpha, ibeta] = df[which(obj == min(obj, na.rm = TRUE))]
      }
    }
  }

  ## Assign to new names
  mat = df.mat
  alpha.array = df.alpha.array
  beta.array = df.beta.array

  ## Recent addition
  if(save){
    cat("Saving aggregated results to ", file.path(destin, resfile), fill=TRUE)
    save(mat,
         alpha.array,
         beta.array,
         df.array,
         file = file.path(destin, resfile))
  }


  ## return(df.mat)
  return(list(mat = mat,
              alpha.array = alpha.array,
              beta.array = beta.array,
              df.array = df.array,
              obj.array = obj.array))
}


##' Helper to aggregate parallelized CV results and obtain the |res| object, all
##' saved in |destin|.
##'
##' @param gridsize Size of CV grid.
##' @param destin Location of saved things.
##'
##' @return List containing, for every (ialpha, ibeta), the "best" result
##'   objects (best in the sense that it had the best likelihood value out of
##'   the |nrep| replicates.)
blockcv_aggregate_res <- function(cv_gridsize, nrep, destin,
                                  sim = FALSE, isim = NULL,
                                  save=FALSE, resfile = "best-res.Rdata"){

  ## df.mat = matrix(NA, ncol=cv_gridsize, nrow=cv_gridsize)
  res.list = list()
  for(ialpha in 1:cv_gridsize){
    for(ibeta in 1:cv_gridsize){

      ## Objective values, over nrep
      obj = rep(NA, nrep)
      ## df = rep(NA, nrep) #
      res.list.inner = list()
      for(irep in 1:nrep){
        ## filename = paste0(ialpha, "-", ibeta, "-", irep, "-fit.Rdata")
        filename = make_refit_filename(ialpha, ibeta, irep, sim, isim)
        tryCatch({
          ## Load fitted result
          load(file.path(destin, filename))
          ## Calculate DF
          res.list.inner[[irep]] = res
          ## Also calculate objective function
          obj[irep] = res$objectives[length(res$objectives)]
        }, error = function(e){})
      }

      ## Calculate the df of the best model
      if(!all(is.na(obj))){
        res.list[[paste0(ialpha, "-", ibeta)]] = res.list.inner[[which.min(obj)]] ## which.min?
      }
    }
  }

  ## Recent addition
  if(save){
    cat("Saving aggregated results to ", file.path(destin, resfile))
    save(res.list, file=file.path(destin, resfile))
  }
  return(res.list)
}

##' An experimental function to visualize the CV results. It calls several other
##' "aggregate" functions in this file.
blockcv_summary <- function(blocktype = 1, datatype = 75, numclust = 5,
                            cv_gridsize = 7,
                            nrep = 5,
                            datadir = "~/Dropbox/research/usc/hpc-output",
                            subfolder = ""){

  ####################
  ## Load data #######
  ####################
  if(is.null(subfolder)) subfolder = ""
  destin = file.path(datadir,
                     paste0("blockcv-", blocktype, "-", datatype, "-", numclust),
                     subfolder)

  ##########################
  ## Get the CV results. ###
  ##########################
  a = blockcv_aggregate(destin, cv_gridsize = cv_gridsize, nfold = 5, nrep = nrep,
                        save=FALSE, resfile = "all-cvres.Rdata")
  cvscore.mat = a$cvscore.mat
  min.inds = a$min.inds

  ## Also get the #nonzero coefficients
  b = blockcv_aggregate_df(gridsize=cv_gridsize, nrep=nrep, destin=destin)
  dfmat = b$mat

  ## Get the refit covarem results
  d = blockcv_aggregate_res(cv_gridsize=cv_gridsize, nrep = nrep, destin=destin)
  res = d[[paste0(min.inds[1] , "-", min.inds[2])]]
  reslist = d
  if(is.null(res)){
    stop(paste0("The model with lambda indices (",
                min.inds[1], ",", min.inds[2], ") is not available."))
  }

  ########################
  ## Get coefficients ####
  ########################
  betalist =  lapply(1:numclust, function(iclust){
    ## Get all betas
    rownames(res$beta[[iclust]])[-1] = colnames(res$X)
    cf = res$beta[[iclust]][-1,, drop=FALSE]
    ## TODO: TRY dplyr here:

    ## Remove the rows that are all zero
    all.zero.rows = which(apply(cf, 1, function(myrow)all(myrow==0)))
    if(length(all.zero.rows) > 0){
      cf = cf[-all.zero.rows,, drop=FALSE]
    }
    round(Matrix::Matrix(cf, sparse=TRUE),3)
  })
  names(betalist) = paste0("Beta matrix, cluster ", 1:numclust)
  pretty.betas = betalist
  colnames(res$alpha)[-1 ] = colnames(res$X)
  alpha = t(res$alpha)
  alpha[which(abs(alpha) < 1E-5)] = 0
  pretty.alphas = round(Matrix::Matrix(alpha, sparse=TRUE),3)

  ######################
  ## Get the sigmas ####
  ######################
  if(res$dimdat == 1){
    pretty.sigmas = sqrt(res$sigma[,1,])
    names(pretty.sigmas) = paste0("Cluster ", 1:numclust)
  } else {
    sigmas = lapply(1:numclust, function(iclust){
      diag(res$sigma[iclust,,])
    })
    names(sigmas) = paste0("Cluster ", 1:numclust)
    pretty.sigmas = lapply(sigmas, sqrt)
  }

  return(list(bestres = res,
              cvscore.mat = cvscore.mat,
              min.inds = min.inds,
              dfmat = dfmat,
              pretty.alphas = pretty.alphas,
              pretty.betas = pretty.betas,
              pretty.sigmas = pretty.sigmas,
              reslist = reslist,
              destin = destin))
}


##' Aggregation wrapper, for simulations
blockcv_summary_sim <- function(nsim = 100,
                                blocktype = 2, datatype = 80, numclust = 2, cv_gridsize = 7,
                                nrep = 5,
                                datadir = "~/Dropbox/research/usc/hpc-output",
                                mc.cores = 1){
  ## datadir = "~/Dropbox/research/usc/hpc-output"
  ## blocktype = 2
  ## datatype = 80
  ## numclust = 2
  ## la('flowcy')
  ## isim = 1
  ## cv_gridsize = 7
  ## nrep = 5
  ## ylist = obj$ylist
  ## obj = generate_data_1d_pseudoreal()
  ## plot_1d(ylist=ylist, res=bestres, countslist=NULL, scale=FALSE)

  ## Form the destin folder
  destin = file.path(datadir,
                     paste0("blockcv-", blocktype, "-", datatype, "-", numclust))


  ## Get |nsim| lists, each containing gridsize^2 best replicates.
  print("Getting all gridsize^2 best replicates, from nsim simulations.")
  start.time = Sys.time()
  reslists = mclapply(1:nsim, function(isim){
    printprogress(isim, nsim, start.time=start.time)
    tryCatch({
      reslist = blockcv_aggregate_res(cv_gridsize = cv_gridsize,
                                      nrep = nrep,
                                      sim = TRUE, isim = isim,
                                      destin = destin)
      return(reslist)
    }, error = function(e){ return(NULL)  })
  }, mc.cores = mc.cores)
  save(reslists, file=file.path(destin, "reslists.Rdata"))
  cat(fill=TRUE)
  print('Saved results to reslist.Rdata')


  ## Get the |min.inds|.
  print("Getting all best CV results, from nsim simulations.")
  start.time = Sys.time()
  cv_info_list = mclapply(1:nsim, function(isim){
    tryCatch({
      printprogress(isim, nsim, start.time=start.time)
      obj = blockcv_aggregate(destin, cv_gridsize, nfold, nrep, sim = TRUE, isim = isim)
      ialpha = obj$min.inds[1]
      ibeta = obj$min.inds[2]
      cvscore = obj$cvscore.mat[ialpha, ibeta]
      return(c(isim = isim, ialpha = ialpha, ibeta = ibeta, cvscore = cvscore))
    }, error=function(e){ return(NULL)  })
  }, mc.cores = mc.cores)
  save(cv_info_list, file=file.path(destin, "cv_info_list.Rdata"))
  cv_info_mat = do.call(rbind, cv_info_list)
  save(cv_info_mat, file=file.path(destin, "cv_info_mat.Rdata"))
  cat(fill = TRUE)
  print('Saved results to cv_info_list.Rdata and cv_info_mat.Rdata')


  ## Get bestres of each of the nsim simulations.
  bestreslist = list()
  for(isim in 1:nsim){
    min.inds = cv_info_mat[isim, c("ialpha", "ibeta")]
    if(!is.null(reslists[[isim]])) next
    reslist = reslists[[isim]]
    bestreslist[[isim]] = reslist[paste0(min.inds[1], "-", min.inds[2])]
  }
  save(bestreslist, file=file.path(destin, "bestreslist.Rdata"))
  print('Saved results to bestreslist.Rdata')


  ## Making a plot of /all/ models
  obj = generate_data_1d_pseudoreal(datadir = "~/stagedir/data")
  ylist = obj$ylist
  print("Making all model plots.")
    start.time = Sys.time()
  reslists = mclapply(1:nsim, function(isim){
    printprogress(isim, nsim, start.time=start.time)
    reslist = reslists[[isim]]
    min.inds = cv_info_mat[isim, c("ialpha", "ibeta")]
    plotname = paste0("sim-", isim, "-", blocktype, "-", datatype, "-", numclust, "-allmodels.png")
    png(file.path(destin, plotname), width = 3000, height = 2000)
    par(mfrow = c(cv_gridsize, cv_gridsize))
    for(ialpha in 1:cv_gridsize){
      for(ibeta in 1:cv_gridsize){
        bestres = reslist[[paste0(ialpha, "-", ibeta)]]
        plot_1d(ylist = ylist, res = bestres,
                countslist = NULL, scale = FALSE, date_axis = FALSE)
        if(all(c(ialpha, ibeta) == min.inds))box(lwd=10,col='blue')
      }
    }
    graphics.off()
    print(paste0("Plot made in ", file.path(destin, plotname)))
  }, mc.cores = mc.cores)
  cat(fill=TRUE)
}


##' Aggregation wrapper, for simulations
blockcv_summary_sim2 <- function(nsim = 100,
                                blocktype = 2, datatype = 80, numclust = 2, cv_gridsize = 7,
                                nrep = 5,
                                outputdir = "~/Dropbox/research/usc/hpc-output",
                                datadir = "~/repos/cruisedat/export",
                                mc.cores = 1, plotonly=FALSE){

  ## Form the destin folder
  destin = file.path(outputdir,
                     paste0("blockcv-", blocktype, "-", datatype, "-", numclust))

  ## If the folder called "summary" doesn't exist, create it.
  if(!dir.exists(file.path(destin, "summary"))){
    dir.create(file.path(destin, "summary"))
  }


  if(!plotonly){
  ## Get |nsim| lists, each containing gridsize^2 best replicates.
  print("Getting all gridsize^2 best replicates, from nsim simulations.")
  start.time = Sys.time()
  reslists = mclapply(1:nsim, function(isim){
    printprogress(isim, nsim, start.time=start.time)
    tryCatch({
      reslist = blockcv_aggregate_res(cv_gridsize = cv_gridsize,
                                      nrep = nrep,
                                      sim = TRUE, isim = isim,
                                      destin = destin)
      return(reslist)
    }, error = function(e){ return(NULL)  })
  }, mc.cores = mc.cores)
  save(reslists, file=file.path(destin, "summary",  "reslists.Rdata"))
  cat(fill=TRUE)
  print('Saved results to reslist.Rdata')



  ## Get the |min.inds|.
  print("Getting all best CV results, from nsim simulations.")
  start.time = Sys.time()
  cv_info_list = mclapply(1:nsim, function(isim){
    tryCatch({
      printprogress(isim, nsim, start.time=start.time)
      obj = blockcv_aggregate(destin, "summary",  cv_gridsize, nfold, nrep, sim = TRUE, isim = isim)
      ialpha = obj$min.inds[1]
      ibeta = obj$min.inds[2]
      cvscore = obj$cvscore.mat[ialpha, ibeta]
      return(c(isim = isim, ialpha = ialpha, ibeta = ibeta, cvscore = cvscore))
    }, error=function(e){ return(NULL)  })
  }, mc.cores = mc.cores)
  save(cv_info_list, file=file.path(destin, "summary",  "cv_info_list.Rdata"))
  cv_info_mat = do.call(rbind, cv_info_list)
  save(cv_info_mat, file=file.path(destin, "summary",  "cv_info_mat.Rdata"))
  cat(fill = TRUE)
  print('Saved results to cv_info_list.Rdata and cv_info_mat.Rdata')


  ## Get bestres of each of the nsim simulations.
  bestreslist = list()
  for(isim in 1:nsim){
    min.inds = cv_info_mat[isim, c("ialpha", "ibeta")]
    if(is.null(reslists[[isim]])) next
    reslist = reslists[[isim]]
    bestreslist[[isim]] = reslist[[paste0(min.inds[1], "-", min.inds[2])]]
  }
  save(bestreslist, file=file.path(destin, "summary",  "bestreslist.Rdata"))
  print('Saved results to bestreslist.Rdata')

  } else {
    ## Load already existing summaries.
    ## load(file=file.path(destin, "summary",  "bestreslist.Rdata"))
    load(file=file.path(destin, "summary",  "reslists.Rdata"))
    load(file=file.path(destin, "summary",  "cv_info_mat.Rdata"))
  }

  ## Making a plot of /all/ models
  if(datatype!=9){
    obj = generate_data_1d_pseudoreal(datadir = datadir)##"~/repos/cruisedat/export")
    ylist = obj$ylist
    countslist = NULL
  } else {
    obj = generate_data_1d_pseudoreal_from_cv(datadir = datadir,##"~/repos/cruisedat/export",
                                              nt1 = 200,
                                              beta_par = 0.3,
                                              p = 10,
                                              bin = TRUE,
                                              dat.gridsize=40)
    ylist = obj$ylist
    countslist = obj$countslist
  }
  print("Making all model plots.")
    start.time = Sys.time()
  mclapply(1:nsim, function(isim){
    printprogress(isim, nsim, start.time=start.time)
    reslist = reslists[[isim]]
    min.inds = cv_info_mat[isim, c("ialpha", "ibeta")]
    plotname = paste0("sim-", isim, "-", blocktype, "-", datatype, "-", numclust, "-allmodels.png")
    png(file.path(destin, "summary",  plotname), width = 3000, height = 2000)
    par(mfrow = c(cv_gridsize, cv_gridsize))
    for(ialpha in 1:cv_gridsize){
      for(ibeta in 1:cv_gridsize){
        bestres = reslist[[paste0(ialpha, "-", ibeta)]]
        ## scale = is.null(countslist)
        plot_1d(ylist = ylist, res = bestres,
                countslist = countslist, scale = scale, date_axis = FALSE)
        if(all(c(ialpha, ibeta) == min.inds))box(lwd=10,col='blue')
      }
    }
    graphics.off()
    print(paste0("Plot made in ", file.path(destin, "summary",  plotname)))
  }, mc.cores = mc.cores)
  cat(fill=TRUE)
}
