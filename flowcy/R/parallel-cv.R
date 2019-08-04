## Synopsis: contains the parallelized main CV wrapper, and all helper functions
## related to cross-validations.

##' CV wrapper for covarem().
##' @param nsplit Number of CV splits. Defaults to 5.
##' @param ... default arguments to covarem().
##' @return List containing (1) the set of coefficients
parallel_cv.covarem <- function(ylist = ylist, X = X,
                                mean_lambdas = NULL,
                                pie_lambdas = NULL,
                                max_mean_lambda = NULL,
                                max_pie_lambda = NULL,
                                gridsize = 9,
                                nsplit = 5,
                                numfork = 3,
                                verbose = FALSE,
                                refit = FALSE,
                                destin = "~",
                                multicore.cv = FALSE,
                                ...){

  ## Printing some information about the parallelism
  if(verbose==TRUE){
    cat("At most ", numfork * nsplit, " cores will be used.", fill = TRUE)
    cat("Output saved to ", destin, fill = TRUE)
  }
  ## TODO: Sangwon: make it so that the warm starts are done by /any/ adjacent
  ## existing point in the grid i.e. code in a search of (+-1, +-1) in the grid.

  ## Basic checks
  stopifnot(length(mean_lambdas) == length(pie_lambdas))
  assert_that(!is.null(max_mean_lambda) | !is.null(mean_lambdas) )
  assert_that(!is.null(max_pie_lambda) | !is.null(pie_lambdas) )
  if(is.null(mean_lambdas)){
    mean_lambdas = c(exp(seq(from = -8, to = log(max_mean_lambda), length = gridsize)))
  }
  if(is.null(pie_lambdas)){
    pie_lambdas = c(exp(seq(from = -8, to = log(max_pie_lambda), length = gridsize)))
  }

  ## Create CV split indices
  assert_that(nsplit >= 2)
  mysplits = cvsplit(ylist, nsplit = nsplit)

  ## Run the CV
  reslist = list()
  start.time = Sys.time()

  ## Define clumps of row numbers (Rows are alpha, columns are beta.)
  ialpha.clumps =
    Map(function(a,b)a:b,
        pmin(seq(from = numfork, to = gridsize+numfork-1, by = numfork), gridsize),
        seq(from = 1, to = gridsize, by = numfork))

  ## Do all of the right edge first.
  ibeta = gridsize
  ialphas = gridsize:1
  move_to_up(ialphas, ibeta,
             pie_lambdas, mean_lambdas,
             gridsize,
             NULL, destin, ylist, X, mysplits, nsplit, refit,
             multicore.cv = multicore.cv,
             ...)

  for(iclump in length(ialpha.clumps):1){

    ialphas = ialpha.clumps[[iclump]]
    cat("clump", iclump, "consists of rows:", ialphas, fill=TRUE)

    ## Traverse from right->left, from the right edge
    new.reslists = mclapply(ialphas, function(ialpha){
      warmstart = loadres(ialpha, gridsize, destin)
      ## cat("Warmstart from (", ialpha, gridsize, ")", fill = TRUE)
      move_to_left(ialpha, (gridsize-1):1,
                   pie_lambdas, mean_lambdas,
                   gridsize,
                   warmstart, destin, ylist, X, mysplits,
                   nsplit, refit,
                   multicore.cv = multicore.cv,
                   ...)
    }, mc.cores = numfork)

    ## ## Start of experimental
    ## onealpha <- function(ialpha){
    ##   warmstart = loadres(ialpha, gridsize, destin)
    ##   ## cat("Warmstart from (", ialpha, gridsize, ")", fill = TRUE)
    ##   move_to_left(ialpha, (gridsize-1):1,
    ##                pie_lambdas, mean_lambdas,
    ##                gridsize,
    ##                warmstart, destin, ylist, X, mysplits,
    ##                nsplit, refit,
    ##                multicore.cv = multicore.cv,
    ##                ...)
    ## }
    ## new.reslists = slurm_apply(ialphas, onealpha,
    ##                            nodes = 2, cpus_per_node = floor(numfork / 2) )
    ## ## End of experimental.

    cat(fill=TRUE)
  }
}


## Should be the same as move_to_left() but with different direction.
move_to_up <- function(ialphas, ibeta,
                       alpha_lambdas, beta_lambdas,
                       gridsize,
                       warmstart, destin,
                       ylist, X, splits, nsplit, refit,
                       multicore.cv = FALSE,
                       ...){
    beginning = TRUE
    assert_that(ibeta == gridsize)
    assert_that(all(diff(ialphas) < 0)) ## Check descending order
    for(ialpha in ialphas){
      cat("(", ialpha, ibeta, ")")
      mywarmstart = (if(beginning){warmstart} else {
                                              ## cat("Warmstart from (", ialpha + 1, ibeta, ")", fill = TRUE)
                                              loadres(ialpha+1, ibeta, destin)
                                            })

      ## Change to cvres!!
      cvres = get_cv_score(ylist, X, splits, nsplit, refit,
                         ## Additional arguments for covarem
                         mean_lambda = beta_lambdas[ibeta],
                         pie_lambda = alpha_lambdas[ialpha],
                         mn = mywarmstart$mn,
                         multicore.cv = multicore.cv,
                         ...)

      ## Get the fitted results on the entire data
      res = covarem(ylist = ylist, X = X,
                    mean_lambda = beta_lambdas[ibeta],
                    pie_lambda = alpha_lambdas[ialpha],
                    mn = mywarmstart$mn, ...)

      saveres(res = res,
              cvres = cvres,
              ialpha = ialpha, ibeta = ibeta, destin = destin,
              beta_lambdas = beta_lambdas,
              alpha_lambdas = alpha_lambdas)
      beginning = FALSE
    }
    cat(fill=TRUE)
}


## Move to the left in a row (fix ialpha)
move_to_left <- function(ialpha, ibetas,
                         alpha_lambdas, beta_lambdas,
                         gridsize,
                         warmstart, destin,
                         ylist, X, splits, nsplit, refit,
                         multicore.cv = FALSE,
                         ...){
  beginning = TRUE
  stopifnot(all(diff(ibetas) < 0)) ## Check descending order
  for(ibeta in ibetas){
    cat("(", ialpha, ibeta, ")")
    mywarmstart = (if(beginning){warmstart} else {
                                            ## cat("Warmstart from (", ialpha, ibeta+1, ")", fill = TRUE)
                                            loadres(ialpha, ibeta+1, destin)
                                          })
      ## mywarmstart =NULL
      cvres = get_cv_score(ylist, X, splits, nsplit, refit,
                           ## Additional arguments for covarem
                           mean_lambda=beta_lambdas[ibeta],
                           pie_lambda=alpha_lambdas[ialpha],
                           mn=mywarmstart$mn,
                           multicore.cv=multicore.cv,
                           ...)

      ## Get the fitted results on the entire data
      res = covarem(ylist=ylist, X=X,
                    mean_lambda=beta_lambdas[ibeta],
                    pie_lambda=alpha_lambdas[ialpha],
                    mn=warmstart$mn, ...)

    saveres(res = res, cvres = cvres, ialpha = ialpha, ibeta = ibeta, destin = destin,
            alpha_lambdas = alpha_lambdas,
            beta_lambdas = beta_lambdas)
    beginning = FALSE
  }
}


## ialpha = 3
## ibeta = 9
##   filename = paste0(ibeta, "-", ialpha, ".Rdata")
##   load(file=file.path(destin, filename))

##' Helper to load parallelized CV results, saved in |destin|.
loadres <- function(ialpha, ibeta, destin){
  filename = paste0(ialpha, "-", ibeta, ".Rdata")
  load(file=file.path(destin, filename))
  return(res)
}

##' Helper to save parallelized CV results, saved in |destin|.
saveres <- function(res, cvres, ialpha, ibeta, destin, alpha_lambdas, beta_lambdas){
  filename = paste0(ialpha, "-", ibeta, ".Rdata")
  save(res, cvres, alpha_lambdas, beta_lambdas, file=file.path(destin, filename))
}


##' Helper to AGGREGATE parallelized CV results, saved in |destin|.
aggregateres <- function(gridsize, destin){

  cvscoremat = matrix(NA, gridsize, gridsize)
  for(ialpha in 1:gridsize){
    for(ibeta in 1:gridsize){
      filename = paste0(ialpha, "-", ibeta, ".Rdata")
      ## Check that results exist
      ## cat("(", ialpha, ibeta, ")", fill=TRUE)
      assert_that(file.exists(file = file.path(destin, filename)))
      ## Load CV score and insert in matrix
      load(file = file.path(destin, filename))
      cvscoremat[ialpha, ibeta] = cvres$mean
    }
  }
  return(cvscoremat)
}


## ## Temporary aggregation
## aggregateres_temp <- function(gridsize, destin, numfork=3){

##   cvscoremat = matrix(NA, gridsize, gridsize)

##   ## Define clumps of row numbers (Rows are alpha, columns are beta.)
##   ialpha.clumps =
##     Map(function(a,b)a:b,
##         pmin(seq(from=numfork, to=gridsize+numfork-1, by=numfork), gridsize),
##         seq(from=1, to=gridsize, by=numfork))

##   for(iclump in 1:length(ialpha.clumps)){

##     ## Fill in right edge first
##     ialphas = ialpha.clumps[[iclump]]
##     ibeta = gridsize
##     for(ialpha in ialphas){
##       filename = paste0(ialpha, "-", ibeta, ".Rdata")
##       cat("(", ialpha, ibeta, ")", fill=TRUE)
##       load(file = file.path(destin, filename))
##       cvscoremat[ialpha, ibeta] = res$mean
##     }

##     ## Traverse from right->left, from the right edge
##     for(ialpha in ialphas){
##       for(ibeta in (gridsize-1):1){
##         cat("(", ialpha, ibeta, ")", fill=TRUE)
##         filename = paste0(ialpha, "-", ibeta, ".Rdata")
##         load(file = file.path(destin, filename))
##         cvscoremat[ialpha, ibeta] = cvres$mean
##       }
##     }
##   }
##   return(cvscoremat)
## }
