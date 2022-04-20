## All directories
outputdir = "/scratch/sangwonh/output"
datadir = "/scratch/sangwonh/data"
scriptdir = "/scratch/sangwonh/scripts/gradients2" ## the current script's location
library(flowmix)

### Read in some command line arguments
parse_args(args = commandArgs(trailingOnly = TRUE), verbose=TRUE)

## Form destination folder
folder1 = paste0("subsample-random-time")
folder2 = paste0("subsample", "-b-", subsample_size)
folder3 = paste0("sim-", isim)
destin = file.path(outputdir, folder1, folder2, folder3)
create_destin(destin)

## Basic setup
maxdev = 0.115525
numclust = 5
cv_gridsize = 7
nrep = 10
nfold = 5

## Load 1d data
dat = readRDS(file = file.path(datadir, "MGL1704-hourly-paper-1d-diam.RDS"))
ylist = dat$ylist
countslist = dat$countslist
X = dat$X
TT = length(ylist)

## CV block size is 10.
blocksize = 20 ## Adjusting for the blocksize

## Use array size to determine the data generation (freeze this!)
nsim0 = nsim
subsamp_ind_list <- lapply(1:nsim0, function(isim0){
  set.seed(isim0 * 10000)
  inds = sample(1:TT, subsample_size, replace = FALSE) %>% sort()
  return(inds)
})
set.seed(NULL)

if(summ){
  cv_summary(destin = destin,
             nfold = nfold,
             nrep = nrep,
             save = TRUE,
             filename = "summary.RDS")
  ## Also save in "subsample-summaries" under destin
  obj = readRDS(file = file.path(destin, "summary.RDS"))
  create_destin(file.path(outputdir, folder1, folder2, "subsample-summaries"))
  q()
}


## Subsample the binned data
inds = subsamp_ind_list[[isim]]
orig_TT = length(ylist)
ylist = ylist[inds]
countslist = countslist[inds]
X = X[inds,]


## Also make the folds
folds = make_cv_folds_subsample_with_original_membership(nfold, blocksize, orig_TT, inds)

## 3. Obtain the maximum regularization parameters
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
  q()
}
prob_lambdas =  logspace(min = 0.0001, max=maxres$alpha, length=cv_gridsize)
mean_lambdas = logspace(min = 0.0001, max=maxres$beta, length=cv_gridsize)


## 4. Estimate model
for(refit in c(FALSE, TRUE)){

  ## Make CV job index matrices
  if(!refit){
    iimat = flowmix::make_iimat(cv_gridsize, nfold, nrep)
  } else {
    iimat = flowmix::make_iimat_small(cv_gridsize)
  }

  ## Divide them up into |arraynum_max| parts
  iilist = make_iilist(arraynum_max, iimat) ## TODO: add flowmix:: here
  if(arraynum > length(iilist)) return(NULL)
  iimat = iimat[iilist[[arraynum]],]

  ## Run the jobs in parallel, in one machine.
  cv.flowmix(
      ## Data
      ylist,
      countslist,
      X,
      ## Define the locations to save the CV.
      destin = destin,
      ## Regularization parameter values
      mean_lambdas,
      prob_lambdas,
      iimat,
      ## Other settings
      maxdev,
      numclust,
      nfold,
      nrep,
      refit = refit,
      save_meta = (arraynum == 1),
      mc.cores = 1,
      blocksize = blocksize,
      folds = folds)
}

