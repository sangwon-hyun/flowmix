####################
## Common setup ####
####################

## Change the output directory to a location of your choice.
outputdir = "~"
datadir = "../paper-data"
library(flowmix)
library(parallel)
mc.cores = 1

## Other setup
maxdev = 0.5
blocksize = 20 ## Size for hourlong blocks (block type 2).
fitonly = FALSE ## Only do refitting?
test = FALSE
sigma_fac = 1

#######################################
## Run one of the next four blocks ####
#######################################


## Figure 8-9: 3d data for a 10-cluster analysis
experiment_type = 1
sim = FALSE
nsim = 1
destin = file.path(outputdir, "3d-10clust")
flowmix::create_destin(destin)
numclust = 10
source("gradients2.R")


## Figure 7: 1d data 5-cluster analysis
experiment_type = 2
sim = FALSE
nsim = 1
numclust = 5
cv_gridsize = 2
nrep = 1
nfold = 2
destin = file.path(outputdir, "1d-5clust")
flowmix::create_destin(destin)
source("gradients2.R")


## Figure 6: Pseudo-real 1d data 5-cluster analysis
experiment_type = 3
destin = file.path(outputdir, "pseudoreal-numclust-simulation")
flowmix::create_destin(destin)
source("gradients2.R")


## For Figure 4-5 (experiment type 4) using artificial 2-cluster 1d simulation data
experiment_type = 4
nsim = 100
for(noise_ii in 1:9){
  destin = file.path(outputdir, "covariate-noise-simulation",
                     paste0("noise-level-", noise_ii))
  flowmix::create_destin(destin)
  for(isim in 1:nsim){
    source("gradients2.R")
  }
}
