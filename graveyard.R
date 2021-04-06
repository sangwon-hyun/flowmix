


##' Aggregation wrapper, for simulations.
##'
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
      print_progress(isim, nsim, start.time=start.time)
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
        print_progress(isim, nsim, start.time=start.time)
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
    print_progress(isim, nsim, start.time=start.time)
    reslist = reslists[[isim]]
    min.inds = cv_info_mat[isim, c("ialpha", "ibeta")]
    plotname = paste0("sim-", isim, "-", blocktype, "-", datatype, "-", numclust, "-allmodels.png")
    grDevices::png(file.path(destin, "summary",  plotname), width = 3000, height = 2000)
    graphics::par(mfrow = c(cv_gridsize, cv_gridsize))
    for(ialpha in 1:cv_gridsize){
      for(ibeta in 1:cv_gridsize){
        bestres = reslist[[paste0(ialpha, "-", ibeta)]]
        scale = !is.null(countslist)
        plot_1d(ylist = ylist, res = bestres,
                countslist = countslist, scale = scale)
        if(all(c(ialpha, ibeta) == min.inds)) graphics::box(lwd = 10, col = 'blue')
      }
    }
    grDevices::graphics.off()
    print(paste0("Plot made in ", file.path(destin, "summary",  plotname)))
  }, mc.cores = mc.cores)
  cat(fill=TRUE)
  }
}




##' Helper to aggregate parallelized CV results and obtain degrees of freedom
##' (DF) estimate, saved in |destin|.
##'
##' @inheritParams
##'
##' @return Matrix containing estimated degrees of freedom.
cv_aggregate_df <- function(destin){


  load(file.path(destin, "meta.Rdata"))
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



##' Temporary plotter.
##'
plot_ylist_other <- function(ylist, countslist, maxcount, col,
                             main = "", cex = 1,
                             ylim = NULL,
                             add = FALSE,
                             bg = NULL){

  stopifnot(ncol(ylist[[1]]) == 1)
  TT = length(ylist)
  if(is.null(ylim)) ylim = range(unlist(ylist))
  if(!add){
    if(!is.null(bg)) par(bg = bg)
    matplot(NA,
            type = 'l',
            lty = 1,
            lwd = .1,
            ylim = ylim,
            xlim = c(1, TT),
            ylab = "",
            xlab = "",
            axes = FALSE)

    ## Add title
    title(main = main)
  }

  ## Handle when countslist is NULL
  ## assertthat::assert_that(!is.null(countslist))
  mx = max(unlist(countslist))
  if(is.null(countslist)) countslist = lapply(ylist, function(y) rep(1, nrow(y)))
  cols = lapply(countslist, function(counts){
    ct = counts / maxcount
    return(sapply(ct, function(one_ct) col %>% adjustcolor(alpha.f = one_ct)))
    ## return(rgb(0, 0, 0, ct))
  })
  pch = 15

  ## Visualize the data
  for(tt in 1:TT){
    y = ylist[[tt]]
    col = cols[[tt]]
    points(x = rep(tt, length(y)),
           y = y, col = col,
           pch = pch, cex = cex)
  }
}



##' A wrapper, to get objective value at a particular time or vector of times
##' \code{times}.
##'
##' @param times Time points of interest (single integer value, or integer vector).
##' @param ... Rest of arguments to \code{objective()}.
##'
##' @return Objective value
objective_subset <- function(times, ...){


  ## Subset things to pass to objective()
  args = list(...)
  args$mu = (args$mu)[times,,,drop=FALSE]
  args$prob = (args$prob)[times,,drop=FALSE]
  args$ylist = (args$ylist)[times]
  args$countslist = (args$countslist)[times]
  stopifnot(all(sapply(args$ylist, nrow) == sapply(args$countslist, length)))

  ## Call the problem
  return(do.call(objective, args))
}



##' Bootstrapping residuals using a coin-flip assignment of cytogram particles.
##'
##' @param ylist_orig Original data, containing the actual particles.
##' @param res A \code{flowmix} class object.
##' @param countslist Counts.
##' @param verbose Print progress if \code{TRUE}. Defaults to \code{FALSE}.
##'
##' @return \code{ylist_bootstrapped} is a dataset that is exactly the same size
##'   as \code{ylist_orig}, but with the particles created by bootstrapping
##'   residuals (after having pooled all of them, across all time points).
##'
##' @examples
##' \dontrun{
##' ## Generate data and fit model
##' set.seed(0)
##' dat1 = generate_data_generic(dimdat = 2, prob1 = 1/8, nt=2000)
##' ylist = dat1$ylist
##' X = dat1$X
##' res = flowmix(ylist, X, numclust=4, niter=300,
##'               mean_lambda = 5E-3, prob_lambda = 5E-3)
##'
##' ## Now, create the bootstrapped datasets.
##' ylists <- lapply(1:nsim, function(isim){
##'   new_ylist = bootstrap(ylist, res)$ylist
##' })
##'
##' ## Make a list of bootstrapped models.
##' bootreslist <- lapply(1:nboot, function(iboot){
##'   ylist = ylists[[iboot]]
##'   res = flowmix(ylist, X, numclust=4, niter=300,
##'                 mean_lambda = 5E-3, prob_lambda = 5E-3)
##' })
##' }
##'
##'
##' @export
##'
bootstrap <- function(ylist, res, countslist = NULL, verbose = FALSE){

  ## Setup
  numclust = res$numclust
  TT = length(ylist)
  dimdat = ylist %>% .[[1]] %>% ncol()

  ## Conduct the E-step once to calculate responsibilities
  resp <- Estep(res$mn, res$sigma, res$prob, ylist = ylist,
                numclust = res$numclust, first_iter = TRUE)

  ## Do a coin-flip draw of membership
  drawslist = draw_membership(resp)
  rm(resp)

  ## Draw residuals
  residuals = get_residuals(ylist, res, drawslist)

  ## Pool all the residuals
  residuals_pooled = lapply(1:res$numclust, function(iclust){
    do.call(rbind, lapply(1:TT, function(tt) residuals[[tt]][[iclust]]))
  })

  ## NEW: Obtain and pool all the biomasses (UGHH)
  if(!is.null(countslist)){

    ## Equivalent  get_residuals()
    index_by_clust <- get_index_by_clust(drawslist)

    counts_pooled = lapply(1:res$numclust, function(iclust){
      do.call(c, lapply(1:TT, function(tt){
        counts = countslist[[tt]]
        index = index_by_clust[[tt]][[iclust]]
        return(counts[index])
      }))
    })

    ## Make sure the number of pooled counts are the same as the number of
    ## pooled residuals
    assertthat::assert_that(all(sapply(counts_pooled, length) == sapply(residuals_pooled, nrow)))
  }


  ## Get the number of coin-flipped particles for each cluster
  ntklist = sapply(residuals, function(resids) resids %>% sapply(., nrow)) %>% t()

  ## Different number of clusters
  ## ntlist = ntklist %>% apply(., 1, sum)
  ## ((1/ntlist) * ntklist ) %>% matplot(type='l')
  ## (ntklist ) %>% matplot(type='l')
  ## plot(ntlist, type='l')
  ## ntlist2 = ylist %>% sapply(., nrow)

  ## Continue here!!! I think it's correct but out of memory... I should just separate out..
  new_y_and_counts_list = list()
  ## <- lapply(1:TT, function(tt){
  for(tt in 1:TT){
    if(verbose) print_progress(tt, TT, fill = TRUE)
    new_y_by_clust <- lapply(1:numclust, function(iclust){
      if(verbose) print_progress(iclust, numclust, fill = TRUE)

      ## How many particles in this cluster to pick?
      ntk = ntklist[tt, iclust]
      if(ntk == 0) return(rep(NA,dimdat) %>% rbind() %>% .[-1,])

      ## Draw the new y's and multiplicities (counts)
      bootrows = sample(1:sum(ntklist[,iclust]), ntk, replace = TRUE)
      resampled_resids <- residuals_pooled[[iclust]][bootrows,,drop=FALSE]
      new_y <- resampled_resids %>% sweep(., 2, res$mn[tt,,iclust], "+")

      ## If there are counts, include
      if(!is.null(countslist)){
        resampled_counts <- counts_pooled[[iclust]][bootrows] ## DO SOMETHING ABOUT THIS
        assertthat::assert_that(length(resampled_counts) == nrow(new_y))
        return(cbind(new_y, counts = resampled_counts))
      } else {
        return(new_y)
      }
    })
    ## return(new_y_by_clust)
    new_y_and_counts_list[[tt]] = new_y_by_clust
  }

  ## Return ylist and countslist separately
  if(!is.null(countslist)){
    new_ylist <- new_y_and_counts_list %>% lapply(., function(a){
      X = do.call(rbind, a) %>% as_tibble()
      return(X %>% select(-counts) %>% as.matrix())
    })
    new_counts <- new_y_and_counts_list %>%
      lapply(., function(a){ do.call(rbind, a)[,"counts", drop=TRUE] %>% unname() })
  } else {
    new_ylist <- new_y_and_counts_list %>% lapply(., function(a){ do.call(rbind, a)})
    new_counts = NULL
  }


  ## saveRDS(list(ylist = new_ylist,
  ##              X = res$X,
  ##              ntklist = ntklist,
  ##              drawslist = drawslist,
  ##              residuals = residuals),
  ##         file = file.path("~/Desktop/bootstrap.RDS"))

  return(list(ylist = new_ylist,
              countslist = new_counts,
              X = res$X,
              ntklist = ntklist,
              drawslist = drawslist,
              residuals = residuals))
}



##' Drawing memberships by coinflip using the "popular vote" i.e. pick
##' membership as cluster with highest responsibility..
##'
##' @param resp List of responsibilities.
##'
##' @export
draw_membership_popular_vote <-function(resp){
  TT = length(resp)
  drawslist = lapply(1:TT, function(tt){
    draws = resp[[tt]] %>% apply(., 1, function(p){
      P = rep(0, length(p))
      P[which.max(p)] = 1
      P
    }) %>% t()
  })
  return(drawslist)
}




##' From subsampling bootstrap results (summary files), produce a set of confidence intervals.
##'
##' @param outputdir Contains files named "summary-(isim).RDS".
##' @param origres Model estimated from the entire dataset.
##' @param ylist_particle Original particle-level cytogram data.
##' @param X Accompanying covariate data.
##'
get_simulated_models <- function(nsim, outputdir, origres, ylist_particle, X, iboot=NULL){

  ## Sample settings
  ## ## 1. Load the original, particle-level data.
  ## obj = readRDS("~/repos/cruisedat/export/MGL1704-hourly-paper-1d-diam-not-binned.RDS")
  ## ylist_particle = obj$ylist
  ## X = obj$X
  ## outputdir = "~/Dropbox/research/usc/hpc-output/subsample-b-148/subsample-summaries"
  ## orig_cvres = cv_summary(destin = file.path("~/Dropbox/research/usc/hpc-output/blockcv-2-76-5"))
  ## origres = orig_cvres$bestres
  ## nsim = 50
  ## get_frequency(nsim, outputdir, orig_cvres, ylist_particle, X, beta = 0.5)
  ## End of temporary


  ## Setup
  numclust = origres$numclust

  ## ## For back-compatibility
  ## origres$prob = origres$pie
  ## class(origres) = "flowmix"

  ## ## TODO: Try to recover nsim directly rom this directory
  ## list.files(outputdir) %>% lapply(., grepl, ...)

  ######################
  ###  Frequencies #####
  ######################
  newres_list = list()
  start.time = Sys.time()
  for(isim in 1:nsim){

    print_progress(isim, nsim, start.time = start.time)

    ## Load the data.
    if(!is.null(iboot)){
      resfile = file.path(outputdir, paste0("summary-", iboot, "-", isim, ".RDS"))
    } else {
      resfile = file.path(outputdir, paste0("summary-", isim, ".RDS"))
    }
    if(!file.exists(resfile)) next
    cvres = readRDS(file = resfile)

    ## Respon
    newres = predict(cvres$bestres, newx = X)
    class(newres) = "flowmix"

    ## Reorder the new res.
    newres = newres %>% reorder_kl(origres, ylist_particle, fac = 100, verbose = FALSE)

    ## Return the results
    newres_list[[isim]] = newres
  }
  return(newres_list)
}
