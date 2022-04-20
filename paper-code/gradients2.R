################################
## Generate data and folds ####
################################
if(experiment_type == 1){

  ## Load the 3d data for a 10-cluster analysis (Figure 8-9)
  ## Data itself was created from the cruisedat R package (demo-average.Rmd)
  dat = readRDS(file = file.path(datadir, "MGL1704-hourly-paper.RDS"))
  list2env(dat, globalenv())

} else if (experiment_type == 2){

  ## Load the 1d data for a 5-cluster analysis (Figure 7)
  ## Data itself was created from the cruisedat R package (demo-average.Rmd)
  dat = readRDS(file = file.path(datadir, "MGL1704-hourly-paper-1d-diam.RDS"))
  list2env(dat, globalenv())

} else if (experiment_type == 3){

  ## Generate realistic 5-cluster 1d data (Figure 6)

  ## 1. Load data (1d-cvres.rds was created from experiment 2)
  source("gradients2-helpers.R")
  obj = generate_data_1d_pseudoreal_from_cv(datadir = datadir) ## loads a file called "1d-cvres.rds"
  ylist = obj$ylist
  X = obj$X
  countslist = obj$countslist

  ## 2. Bin the data
  dat.grid = flowmix::make_grid(ylist, gridsize = 40)
  obj = flowmix::bin_many_cytograms(ylist, dat.grid, mc.cores = 1, verbose = TRUE)
  ylist = lapply(obj$ybin_list, cbind)
  countslist = obj$counts_list

  ## Other setup
  maxdev = 0.1155245
  final_range = lapply(ylist, range) %>% unlist %>% range
  sigma_fac = (final_range[2] - final_range[1])/numclust

} else if (experiment_type == 4){

  ## Generate artificial 2-cluster 1d simulation data, in a few steps (Figure 4-5)

  ## 1. Generate fake 1d data
  obj = generate_data_1d_pseudoreal(bin = FALSE, datadir = datadir, nt1 = 200, beta_par = 0.3, p = 10)
  X = obj$X
  ylist = obj$ylist
  countslist = obj$countslist
  TT = length(ylist)

  ## 2. Bin with just counts
  dat.grid = flowcy::make_grid(ylist, gridsize = 40)
  obj = flowcy::bin_many_cytograms(ylist, dat.grid, mc.cores = 1, verbose = TRUE)
  ylist = lapply(obj$ybin_list, cbind)
  countslist = obj$counts_list

  ## Check about noise level
  sigmalist = seq(from = 3 * 0.1, to = 3 * 0.9, length = 9)
  stopifnot(noise_ii %in% c(1:9))
}


################################
## Calculate maximum lambda ####
################################
if(file.exists(file.path(destin, "maxres.Rdata"))){
  load(file.path(destin, "maxres.Rdata"), verbose = TRUE)
} else {
  maxres = get_max_lambda(destin, "maxres.Rdata",
                          ylist = ylist,
                          countslist = countslist,
                          X = X,
                          numclust = numclust,
                          maxdev = maxdev,
                          max_prob_lambda = 1,
                          max_mean_lambda = 3,
                          verbose = TRUE)
}
prob_lambdas =  logspace(min = 0.0001, max=maxres$alpha, length=cv_gridsize)
mean_lambdas = logspace(min = 0.0001, max=maxres$beta, length=cv_gridsize)


## ## Temporary: time the code for small prob_lambda values
## library(profvis)
## library(magrittr)
## library(stringr)
## library(tidyverse)
## la('flowmix')

## mstep_beta_time = total_time = fortran_time = matrix(0, nrow=cv_gridsize, ncol=cv_gridsize)
## for(ialpha in 1:cv_gridsize){
##   for(ibeta in 1:cv_gridsize){
##     outfile = file.path(destin, paste0("profvis-", ialpha, "-", ibeta, ".out"))
##     ## if(!file.exists(outfile)){
##       obj = profvis::profvis({
##         flowmix(ylist=ylist, countslist=countslist,
##                 X = X,
##                 numclust = 5,
##                 prob_lambda = prob_lambdas[ialpha],
##                 mean_lambda = mean_lambdas[ibeta],
##                 niter = 20,
##                 nrep = 1,
##                 verbose = FALSE)
##       }, prof_output = outfile)
##     ## }
##     times <- outfile %>% summaryRprof() %$% by.total %>% select(total.time)
##     total_time[ialpha, ibeta] <- times %>% rownames_to_column("type") %>% filter(str_detect(type, 'profvis')) %>% select(total.time) %>% unlist()
##     fortran_time[ialpha, ibeta] <- times %>% rownames_to_column("type") %>% filter(str_detect(type, '.Fortran')) %>% select(total.time) %>% unlist()
##     mstep_beta_time[ialpha, ibeta] <- times %>% .['\"admm_oneclust\"',] %>% unlist()
##     la_admm_time[ialpha, ibeta] <- times %>% .['\"la_admm_oneclust\"',] %>% unlist()
##   }
## }

## rownames(fortran_time) = rownames(total_time) = rownames(mstep_beta_time) = prob_lambdas %>% signif(3)
## colnames(fortran_time) = colnames(total_time) = rownames(mstep_beta_time) = prob_lambdas %>% signif(3)

## #### Total time spent ###########################################################################################
## total_time %>% drawmat_precise(xlab = expression(lambda[beta]),
##                                ylab = expression(lambda[alpha]),
##                                main="Total time of for flowmix()")
## par(mfrow=c(2,1), mar=c(4,5,1,3))
## matplot(x=prob_lambdas, y=total_time, type='l', lwd=.5, col='red', lty=1, ylab="Total time of flowmix()",
##         pch=16, cex=.5, xlab = expression(lambda[alpha]))
## abline(v=prob_lambdas, col=rgb(0,0,0,0.5), lwd=0.3)
## matplot(x=prob_lambdas, y=total_time, type='l', lwd=.5, col='red', lty=1, ylab="Total time of flowmix()",
##         log = "x", pch=16, cex=.5, xlab = expression(lambda[alpha]))
## abline(v=prob_lambdas, col=rgb(0,0,0,0.5), lwd=0.3)

## #### Time spent by Fortran #####################################################################################
## fortran_time %>% drawmat_precise(xlab = expression(lambda[beta]),
##                                  ylab = expression(lambda[alpha]),
##                                  main="Total time spent by \n .Fortran()")
## (fortran_time / total_time) %>% drawmat_precise(xlab = expression(lambda[beta]),
##                                                 ylab = expression(lambda[alpha]),
##                                                 main="Proportion of time spent by\n .Fortran()")
## par(mfrow=c(3,1), mar=c(4,5,1,3))
## fortran_time %>% matplot(x=prob_lambdas, y=., type='l', lwd=.5, col='red', lty=1, ylab="Total time spent by .Fortran()",
##                          pch=16, cex=.5, xlab = expression(lambda[alpha]))
## abline(v=prob_lambdas, col=rgb(0,0,0,0.5), lwd=0.3)
## fortran_time %>% matplot(x=prob_lambdas, y=., type='l', lwd=.5, col='red', lty=1, ylab="Total time spent by .Fortran()",
##         log = "x", pch=16, cex=.5, xlab = expression(lambda[alpha]))
## abline(v=prob_lambdas, col=rgb(0,0,0,0.5), lwd=0.3)
## (fortran_time/total_time) %>% matplot(x=prob_lambdas, y=., type='l', lwd=.5, col='red', lty=1, ylab="Proportion of time spent by .Fortran()",
##         log = "x", pch=16, cex=.5, xlab = expression(lambda[alpha]), ylim=c(0,1))
## abline(v=prob_lambdas, col=rgb(0,0,0,0.5), lwd=0.3)

## #### Total time spent by beta M step ############################################################################
## mstep_beta_time %>% drawmat_precise(xlab = expression(lambda[beta]),
##                                     ylab = expression(lambda[alpha]),
##                                     main="Total time of beta M step")
## (mstep_beta_time/total_time) %>% drawmat_precise(xlab = expression(lambda[beta]),
##                                     ylab = expression(lambda[alpha]),
##                                     main="Proportion of time of beta M step")

## par(mfrow=c(2,1), mar=c(4,5,1,3))
## t(mstep_beta_time) %>% matplot(x=prob_lambdas, y=., type='l', lwd=.5, col='red', lty=1, ylab="Total time spent by \n Beta M step",
##                                          pch=16, cex=.5, xlab = expression(lambda[beta]), log="x",
##                                          main="total time of beta M step")
## abline(v=prob_lambdas, col=rgb(0,0,0,0.5), lwd=0.3)
## t(mstep_beta_time/total_time) %>% matplot(x=mean_lambdas, y=., type='l', lwd=.5, col='red', lty=1, ylab="Proportion of time spent by \n Beta M step",
##                                          pch=16, cex=.5, xlab = expression(lambda[beta]), log="x",
##                                          main="Proportion of time of beta M step",
##                                          ylim = c(0,1))
## abline(v=prob_lambdas, col=rgb(0,0,0,0.5), lwd=0.3)


######################
## Create CV folds ###
######################
folds = make_cv_folds(ylist, nfold, blocksize = blocksize)

############################
### Save meta data once ####
############################
save(folds,
     nfold,
     nrep, ## Added recently
     cv_gridsize,
     mean_lambdas,
     prob_lambdas,
     ylist, countslist, X,
     ## Save the file
     file = file.path(destin, 'meta.Rdata'))
print(paste0("wrote meta data to ", file.path(destin, 'meta.Rdata')))


###################################################
## Run ONE CV job (ialpha, ibeta, ifold, irep) ####
###################################################
iimat = make_iimat(cv_gridsize, nfold, nrep)

## Run the jobs in parallel, in one machine.
ind = iimat[,"ind"]
parallel::mclapply(ind, function(ii){

  ## Indices defining the jobs to run
  ialpha = iimat[ii,"ialpha"]
  ibeta = iimat[ii,"ibeta"]
  ifold = iimat[ii,"ifold"]
  irep = iimat[ii,"irep"]

  cat('(ialpha, ibeta, ifold, irep)=', c(ialpha, ibeta, ifold, irep), fill = TRUE)

  for(isim in 1:nsim){

    ## (for experiment type 4, Figure 5) Add noise to X
    if(experiment_type == 4){
      sigmalist = seq(from = 3 * 0.1, to = 3 * 0.9, length = 9)
      X = add_noise(X, sigmalist, noise_ii)
    }

    ## Test
    one_job(ialpha = ialpha,
            ibeta = ibeta,
            ifold = ifold,
            irep = irep,
            folds = folds,
            destin = destin,
            mean_lambdas = mean_lambdas,
            prob_lambdas = prob_lambdas,
            ## Arguments for covarem()
            ylist = ylist, countslist = countslist, X = X,
            sigma_fac = sigma_fac,
            ## Additional arguments for covarem(), for ellipsis.
            numclust = numclust,
            maxdev = maxdev,
            ## Simulations? Yes or no
            sim = sim,
            isim = isim,
            verbose = TRUE)
  }
  return(NULL)
}, mc.cores = mc.cores)


###################################################
## Run ONE CV job (ialpha, ibeta, ifold, irep) ####
###################################################
## Define all the jobs.
iimat = make_iimat_small(cv_gridsize)

parallel::mclapply(iilist, function(ii){

  ialpha = iimat[ii, "ialpha"]
  ibeta = iimat[ii, "ibeta"]
  print(paste0('refitting (ialpha, ibeta)=(', ialpha, ",", ibeta, ")"))

  for(isim in 1:nsim){

    ## (for experiment type 4, Figure 5) Add noise to X
    if(experiment_type == 4){
      sigmalist = seq(from = 3 * 0.1, to = 3 * 0.9, length = 9)
      X = add_noise(X, sigmalist, noise_ii)
    }

    one_job_refit(ialpha = ialpha, ibeta = ibeta, destin = destin,
                  mean_lambdas = mean_lambdas, pie_lambdas = pie_lambdas,
                  ## Arguments to covarem()
                  ylist = ylist, countslist = countslist, X = X,
                  sigma_fac = sigma_fac,
                  numclust = numclust,
                  maxdev = maxdev,
                  nrep = nrep,
                  sim = sim,
                  isim = isim)
  }
  return(NULL)
}, mc.cores = mc.cores)


load(file=file.path(destin, "meta.Rdata"), verbose = TRUE)
obj = profvis({
  set.seed(0)
  res = flowmix(numclust = 5, maxdev = 0.5,
                verbose = TRUE, ylist = ylist, countslist = countslist,
                X = X, mean_lambda = 1E-3,
                prob_lambda = 1E-2, ## Try 1E-2, 1E-4, 1E-6, 1E-7
                niter = 20, nrep = 1)
})
print(obj)
