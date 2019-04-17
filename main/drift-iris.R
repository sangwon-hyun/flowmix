## Synopsis: Try code on artificial IRIS data.

## Setup
la("~/repos/flowcy/flowcy")
source("~/repos/flowcy/main/drift-helpers.R")
library(MASS)
library(testthat)
plotdir = "/home/shyun/Desktop/cytometry-plots"

## Load data:
## load("data.Rdata")                      #Placeholder
dimdat = 3
numclust = 2
T = 3
## data = make_data(T)
data = make_data(T, fac=0.5)
## par(mfrow=c(1,3))
## sapply(data, function(a)plot(a[,2:3], col=iris[,5], pch=16, cex=3))

## Initialize clusters
la("~/repos/flowcy/flowcy")
pie = init_pi(numclust, T)
pie = rbind(c(0.2,0.8), c(0.3, 0.7), c(0.7, 0.3))
mu = init_mu(data, numclust, T)
sigma = init_sigma(data, numclust, T, fac=1)
n = nrow(data[[1]])
s = 1E-4 #1E-4
lam1 = 0 ## lam2 = 10


## All outer steps
for(lam2 in c(0,1,2,5,6,7,8,9,10,15,20)){
## for(lam2 in c(6,7,8,9)){
  lam2 = 100

  results = driftem(data, mu, pie, niter=1000, sigma, T, tol1 = 1E-10, tol2 = 1E-4, lam1, lam2, s)

  ## Visualize objective values
  filename = file.path(plotdir,
                paste0("objectives-drifts-lam1-", lam1, "-lam2-", lam2, ".pdf"))
  pdf(filename, width=5, height=5)
  plot(results$objectives, type='l', lwd=2)
  graphics.off()

  ## Visualize the cluster centers
  filename = file.path(plotdir,
                paste0("drifts-new-lam1-", lam1, "-lam2-", lam2, ".pdf"))
  print(filename)
  pdf(filename, width=8, height=8)
  myplot(data, results$mulist[[results$final.iter]], lam1, lam2)
  graphics.off()
}


## Trying to see what happens if smaller step size (s) is taken
filename = file.path(plotdir,
              paste0("drifts-new-lam1-", lam1, "-lam2-", lam2, "-iterations", ".pdf"))
print(filename)
pdf(filename, width=8, height=8)
for (iter in round(seq(from=1, to=results$final.iter, length=min(100, results$final.iter)))){
  myplot(data, results$mulist[[iter]], lam1, lam2)
  title(sub=iter)
}
graphics.off()

  plot(results$objectives, type='l', lwd=2)
## Going over the iters$.


## Trying to see Rprof
la("~/repos/flowcy/flowcy")
lam2 = 10
Rprof(tmp <- tempfile())                # Starts profiling
results = driftem(data, mu, pie, niter=4, sigma, T, tol1 = 1E-10, tol2 = 1E-4, lam1, lam2, s, numclust)
Rprof()                                 # End profiling
b = summaryRprof(tmp)
a
b
## Okay, so there is evidence that the changing of the order of the loop helped.


## Trying to fix pie step
la("~/repos/flowcy/flowcy")
s = 1E-4
lam1 = 1000
lam2 = 0
results = driftem(data, mu, pie, niter=1000, sigma, T, tol1 = 1E-10, tol2 = 1E-4, lam1, lam2, s)
results$sigma
source("~/repos/flowcy/main/drift-helpers.R")
myplot(data, results$mulist[[results$final.iter]], results$pielist[[results$final.iter]],
       results$sigmalist[[results$final.iter]],
       lam1, lam2)
## the pi



## Seeing what happens when drift EM
la("~/repos/flowcy/flowcy")
results = driftem(data, mu, pie, niter=1000, sigma, T,
                  tol1 = 1E-10, tol2 = 1E-4, lam1, lam2, s, numclust)
## Compare the results.
