
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
  obj = flowmix::generate_data_1d_pseudoreal_from_cv(datadir = datadir) ## loads a file called "1d-cvres.rds"
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
     pie_lambdas,
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
