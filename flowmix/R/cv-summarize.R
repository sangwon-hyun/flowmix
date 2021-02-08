##' Main function for summarizing the cross-validation results.
##'
##' @inheritParams cv.flowmix
##' @param filename File name to save to.
##'
##' @return List containing various outcomes from the cross-validation, such as
##'   \code{bestres} which is the \code{flowmix} class object of the overall
##'   best model chosen from the cross-validation; \code{cvscoremat} containing
##'   a 2d matrix of CV scores from all pairs of lambdas; \code{bestreslist}
##'   contains all the best models (out of \code{nrep} EM replications} from the
##'   each pair of lambda values. If \code{isTRUE{save}}, nothing is returned.
##'
##' @export
cv_summary <- function(destin = ".",
                       cv_gridsize,
                       nrep,
                       nfold,
                       save = FALSE,
                       filename = "summary.RDS"
                       ){

  ####################
  ## Load data #######
  ####################
  load(file = file.path(destin, 'meta.Rdata'), verbose = FALSE)

  ##########################
  ## Get the CV results. ###
  ##########################
  a = cv_aggregate(destin, cv_gridsize = cv_gridsize, nfold = nfold, nrep = nrep)
  cvscore.mat = a$cvscore.mat
  min.inds = a$min.inds

  ## Also get the #nonzero coefficients (not used now)
  dfmat = cv_aggregate_df(cv_gridsize = cv_gridsize, nrep = nrep, destin = destin)$mat

  ## Get the refit flowmix results
  bestreslist = cv_aggregate_res(cv_gridsize = cv_gridsize, nrep = nrep, destin = destin)
  bestres = bestreslist[[paste0(min.inds[1] , "-", min.inds[2])]]
  if(is.null(bestres)){
    stop(paste0("The model with lambda indices (",
                min.inds[1], ",", min.inds[2], ") is not available."))
  }

  if(is.null(colnames(bestres$X))){
    colnames(bestres$X) = 1:ncol(bestres$X)
  }

  ########################
  ## Get coefficients ####
  ########################
  betalist =  lapply(1:bestres$numclust, function(iclust){
    ## Get all betas
    rownames(bestres$beta[[iclust]])[-1] = colnames(bestres$X)
    cf = bestres$beta[[iclust]][-1,, drop=FALSE]
    ## TODO: TRY dplyr here:

    ## Remove the rows that are all zero
    all.zero.rows = which(apply(cf, 1, function(myrow)all(myrow == 0)))
    if(length(all.zero.rows) > 0){
      cf = cf[-all.zero.rows,, drop=FALSE]
    }
    round(Matrix::Matrix(cf, sparse=TRUE),3)
  })
  names(betalist) = paste0("Beta matrix, cluster ", 1:bestres$numclust)
  pretty.betas = betalist
  colnames(bestres$alpha)[-1 ] = colnames(bestres$X)
  alpha = t(bestres$alpha)
  alpha[which(abs(alpha) < 1E-5)] = 0
  pretty.alphas = round(Matrix::Matrix(alpha, sparse=TRUE),3)

  ######################
  ## Get the sigmas ####
  ######################
  if(bestres$dimdat == 1){
    pretty.sigmas = sqrt(bestres$sigma[,1,])
    names(pretty.sigmas) = paste0("Cluster ", 1:bestres$numclust)
  } else {
    sigmas = lapply(1:bestres$numclust, function(iclust){
      diag(bestres$sigma[iclust,,])
    })
    names(sigmas) = paste0("Cluster ", 1:bestres$numclust)
    pretty.sigmas = lapply(sigmas, sqrt)
  }

  out = list(bestres = bestres,
             cvscore.mat = cvscore.mat,
             min.inds = min.inds,
             dfmat = dfmat,
             ## Pretty formatted data
             pretty.alphas = pretty.alphas,
             pretty.betas = pretty.betas,
             pretty.sigmas = pretty.sigmas,
             ## List of all best models for all lambda pairs.
             bestreslist = bestreslist,
             destin = destin)

  if(save){ saveRDS(out, file=file.path(destin, filename)); return(NULL) }
  return(out)
}

##' Aggregate CV scores from the results, saved in \code{destin}.
##'
##' @inheritParams cv.flowmix
##'
##' @export
cv_aggregate <- function(destin, cv_gridsize, nfold, nrep,
                         sim = FALSE, isim = 1){

  ## ## Read the meta data (for |nfold|, |cv_gridsize|, |nrep|)
  load(file = file.path(destin, 'meta.Rdata'))

  ## For back-compatibility
  if(exists("pie_lambda")) prob_lambda = pie_lambda
  if(exists("pie_lambdas")) prob_lambdas = pie_lambdas

  ## Aggregate the results
  cvscore.array = array(NA, dim = c(cv_gridsize, cv_gridsize, nfold, nrep))
  cvscore.mat = matrix(NA, nrow = cv_gridsize, ncol = cv_gridsize)
  for(ialpha in 1:cv_gridsize){
    for(ibeta in 1:cv_gridsize){
      obj = matrix(NA, nrow=nfold, ncol=nrep)
      for(ifold in 1:nfold){
        for(irep in 1:nrep){
          filename = make_cvscore_filename(ialpha, ibeta, ifold, irep, sim, isim)
          tryCatch({
            load(file.path(destin, filename), verbose = FALSE)

            ## Purely for back-compatability
            if(exists("pie_lambda")) prob_lambda = pie_lambda
            if(exists("pie_lambdas")) prob_lambdas = pie_lambdas

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
  rownames(cvscore.mat) = signif(prob_lambdas,3)
  colnames(cvscore.mat) = signif(mean_lambdas,3)


  ## Find the minimum
  mat = cvscore.mat
  min.inds = which(mat == min(mat, na.rm=TRUE), arr.ind=TRUE)

  ## Return the results
  out = list(cvscore.array = cvscore.array,
              cvscore.mat = cvscore.mat,
              mean_lambdas = mean_lambdas,
              prob_lambdas = prob_lambdas,
              min.inds = min.inds)
  return(out)
}




##' Helper to aggregate parallelized CV results and obtain the |res| object, all
##' saved in |destin|.
##'
##' @inheritParams cv_gridsize
##' @inheritParams nrep
##' @inheritParams destin
##'
##' @return List containing, for every (ialpha, ibeta), the "best" estimated
##'   model out of the |nrep| replicates (best in the sense that it had the best
##'   likelihood value out of the |nrep| replicates.)
cv_aggregate_res <- function(cv_gridsize, nrep, destin
                             ## Is this a simulation or not? (soon to be outdated)
                             ## sim = FALSE,
                             ## isim = NULL
                             ){

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
        filename = make_refit_filename(ialpha, ibeta, irep)##, sim, isim)
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
  return(res.list)
}



##' Aggregation wrapper, for simulations.
##'
##' @inheritParams cv_gridsize
##' @inheritParams nfold
##' @inheritParams nrep
##' @param blocktype block type
##' @param datatype data type
##' @param numclust number of clusters
##' @param outputdir output directory, where the destination folder exists.
##'
##' @export
cv_summary_sim <- function(nsim = 100,
                           blocktype = 2, datatype = 80, numclust = 2, cv_gridsize = 7,
                           nrep = 5,
                           outputdir = "~/Dropbox/research/usc/hpc-output",
                           datadir = "~/repos/cruisedat/export",
                           mc.cores = 1,
                           plotonly = FALSE,
                           plot_all_models = FALSE){

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
        reslist = cv_aggregate_res(cv_gridsize = cv_gridsize,
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
        obj = cv_aggregate(destin, cv_gridsize = cv_gridsize, nfold = nfold, nrep = nrep, sim = TRUE, isim = isim)##, save=FALSE)
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
      print(bestreslist)
    print('Saved results to bestreslist.Rdata')

  } else {
    ## Load already existing summaries.
    load(file=file.path(destin, "summary",  "bestreslist.Rdata"))
    load(file=file.path(destin, "summary",  "reslists.Rdata"))
    load(file=file.path(destin, "summary",  "cv_info_mat.Rdata"))
  }

  ## Making a plot of /all/ models
  if(plot_all_models){
  if(datatype!=9){
    obj = generate_data_1d_pseudoreal(datadir = datadir,
                                      nt = 200,
                                      beta_par = 0.3,
                                      p = 10,
                                      bin = TRUE,
                                      dat.gridsize = 40)##"~/repos/cruisedat/export")
    ylist = obj$ylist
    countslist = obj$countslist
  } else {
    obj = generate_data_1d_pseudoreal_from_cv(datadir = datadir)##"~/repos/cruisedat/export",
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
        scale = !is.null(countslist)
        plot_1d(ylist = ylist, res = bestres,
                countslist = countslist, scale = scale, date_axis = FALSE)
        if(all(c(ialpha, ibeta) == min.inds))box(lwd=10,col='blue')
      }
    }
    grDevices::graphics.off()
    print(paste0("Plot made in ", file.path(destin, "summary",  plotname)))
  }, mc.cores = mc.cores)
  cat(fill=TRUE)
  }
}



##' From all the folds/replicates, get rid of the |nrep| replicates by retaining
##' only the best one. (NOT USED NOW)
cv_reduce_by_nrep <- function(destin, cv_gridsize, nfold, nrep, sim = FALSE,
                              isim = NULL, save = FALSE,
                              resfile = "all-cvres.Rdata"){

  ## Read the meta data (for |nfold|, |cv_gridsize|, |nrep|)
  load(file = file.path(destin, 'meta.Rdata'))

  ## Identify the best CV replicate
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
  rownames(cvscore.mat) = signif(prob_lambdas,3)
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
         prob_lambdas,
         min.inds,
         file = file.path(destin, resfile))
  }

  return(list(cvscore.array = cvscore.array,
              cvscore.mat = cvscore.mat,
              mean_lambdas = mean_lambdas,
              prob_lambdas = prob_lambdas,
              min.inds = min.inds))
}



##' Helper to aggregate parallelized CV results and obtain degrees of freedom
##' (DF) estimate, saved in |destin|.
##'
##' @inheritParams gridsize
##' @inheritParams destin
##' @inheritParams nrep
##'
##' @return Matrix containing estimated degrees of freedom.
cv_aggregate_df <- function(cv_gridsize, nrep, destin){

  df.array = obj.array = df.alpha.array = df.beta.array = array(NA, dim=c(cv_gridsize, cv_gridsize, nrep))
  df.mat = df.alpha.mat = df.beta.mat = matrix(NA, ncol=cv_gridsize, nrow=cv_gridsize)
  for(ialpha in 1:cv_gridsize){
    for(ibeta in 1:cv_gridsize){

      ## Objective value
      obj = rep(NA, nrep)
      df = df.alpha = df.beta = rep(NA, nrep)
      for(irep in 1:nrep){

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
        min.df = df[which(obj == min(obj, na.rm = TRUE))]
        df.mat[ialpha, ibeta] = min.df[1] ## RARELY there are duplicates.. (especially when bootstrap is done)
      }
    }
  }

  ## Assign to new names
  mat = df.mat
  alpha.array = df.alpha.array
  beta.array = df.beta.array

  ## return(df.mat)
  out = list(mat = mat,
             alpha.array = alpha.array,
             beta.array = beta.array,
             df.array = df.array,
             obj.array = obj.array)
  return(out)
}
