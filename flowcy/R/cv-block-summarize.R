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


##' Apply the 1SE rule (results are similar to those of get_optimal_info())
##' @param outputdir Location of the output files. e.g. outputdir="~/output/vanilla"
##' @return NULL.
blockcv_onese <- function(destin,
                          gridsize,
                          nrep,
                          sim=FALSE, isim=NULL){
  ## outputdir="~/repos/flowcy/output/vanilla"

  ## Gather the degrees of freedom
  ## dfmat = aggregateres_df(gridsize, destin)
  ## cvscoremat = aggregateres(gridsize, destin)
  cvscore.mat = blockcv_aggregate(destin = destin,
                                 cv_gridsize = cv_gridsize,
                                 nfold = nfold,
                                 sim = TRUE,
                                 isim = 1,
                                 nrep = nrep)$cvscore.mat

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
    cat("Saving aggregated results to ", file.path(destin, resfile))
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
                                datadir = "~/Dropbox/research/usc/hpc-output"){
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
  reslists = lapply(1:nsim, function(isim){
    printprogress(isim, nsim)
    reslist = blockcv_aggregate_res(cv_gridsize = cv_gridsize, nrep = nrep,
                                    sim = TRUE, isim = isim,
                                    destin = destin)
    ## min.inds = reslist
    ## bestres = reslist[[paste0(min.inds[1], "-", min.inds[2])]]
    ## return(bestres)
    return(reslist)
  })

  ## Making a plot of /all/ models
  png(file.path(destin, paste0(blocktype, "-", datatype, "-", numclust, "-allmodels.png")),
      width = 3000, height = 2000)
  par(mfrow=c(cv_gridsize, cv_gridsize))
  for(ialpha in 1:cv_gridsize){
    for(ibeta in 1:cv_gridsize){
      bestres = reslist[[paste0(ialpha, "-", ibeta)]]
      plot_1d(ylist = ylist, res = bestres,
              countslist = NULL, scale = FALSE, date_axis = FALSE)
    }
  }
  graphics.off()

  ## Get the |min.inds|.
  print("Getting all best CV results, from nsim simulations.")
  cv_info_list = lapply(1:nsim, function(isim){
    printprogress(isim, nsim)
    obj = blockcv_aggregate(destin, cv_gridsize, nfold, nrep, sim = TRUE, isim = isim)
    ialpha = obj$min.inds[1]
    ibeta = obj$min.inds[2]
    cvscore = obj$cvscore.mat[ialpha, ibeta]
    c(isim = isim, ialpha = ialpha, ibeta = ibeta, cvscore = cvscore)
  })
  cv_info_mat = do.call(rbind, cv_info_list)

  ## Get each of the best guys.
  bestreslist = list()
  for(isim in 1:nsim){
    min.inds = min.inds.mat[isim, c("ialpha", "ibeta")]
    reslist = reslists[[isim]]
    bestreslist[[isim]] = reslist[paste0(min.inds[1], "-", min.inds[2])]
  }
  return(bestreslist)
  ## return(list(bestreslist = bestreslist))
}
