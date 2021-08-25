####################
## Common setup ####
####################


## Change the output directory to a location of your choice.
outputdir = "~"
datadir = "../paper-data"

library(flowmix) ## Caution! use the version of the package as of Mar 2021.
library(parallel)
mc.cores = 1

## Source in some helper functions

## Other setup
maxdev = 0.5
blocksize = 20 ## Size for hourlong blocks (block type 2).
test = FALSE
sigma_fac = 1

#######################################
## Run one of the next four blocks ####
#######################################


## Figure 8-9: 3d data for a 10-cluster analysis
experiment_type = 1
sim = FALSE
nsim = 1
numclust = 10
cv_gridsize = 10
nrep = 10
nfold = 5
destin = file.path(outputdir, "3d-10clust")
flowmix::create_destin(destin)
## maxres = list(alpha = 0.125, beta = 0.375)
source("gradients2.R")


## Figure 7: 1d data 5-cluster analysis
## (Figure 7, 12, 13, and Table 2 are from here).
experiment_type = 2
sim = FALSE
nsim = 1
numclust = 5
cv_gridsize = 7
nrep = 10
nfold = 5
maxdev = 0.115525
destin = file.path(outputdir, "1d-5clust")
flowmix::create_destin(destin)
## maxres = list(alpha= 0.5, beta = 1.5)
source("gradients2.R")


## Figure 6: Simulation using pseudo-real 1d data 5-cluster analysis
experiment_type = 3
sim = TRUE
nsim = 100
for(numclust in 2:8){
  destin = file.path(outputdir, "pseudoreal-numclust-simulation",
                     paste0("numclust-", numclust))
  flowmix::create_destin(destin)
  source("gradients2.R")
}


## For Figure 4-5 Simulation using artificial 2-cluster 1d synthetic data
experiment_type = 4
sim = TRUE
nsim = 100
for(noise_ii in 1:9){
  destin = file.path(outputdir, "covariate-noise-simulation",
                     paste0("noise-level-", noise_ii))
  flowmix::create_destin(destin)
  for(isim in 1:nsim){
    source("gradients2.R")
  }
}
