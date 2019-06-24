context("Test optimization for ths drifting clusters problem")

test_that("Initialization works well.", {
  data(iris)
  dimdat = 3
  datapoints = iris[,1:dimdat]
  numclust = 2
  T = 3
  pie = init_pi(numclust, T)
  mu = init_mu(data, numclust, T)
  sigma = init_sigma(data, numclust, T)
  testthat::expect_equal(dim(pie), c(T, numclust))
  testthat::expect_equal(dim(mu), c(T, numclust, dimdat))
  testthat::expect_equal(dim(sigma), c(T, numclust, dimdat, dimdat))
})


test_that("Mstep_mu_exact() works well.", {
Mstep_mu_exact

  testthat::expect_equal(dim(sigma), c(T, numclust, dimdat, dimdat))


    ## ## Makepmat needs to be
    ## Pmat = makePmat(dimdat, TT)
})


test_that("EM works well.", {

  pielist = mulist = sigmalist = list()
  objectives = c()
  pielist[[1]] = pie
  mulist[[1]] = mu
  objectives[[1]] = objective_overall(mulist[[1]], pielist[[1]], sigma, data)
  if(class(data[[1]])!="matrix"){  data = lapply(data, as.matrix) }

  ## Load data:
  ## load("data.Rdata")                      #Placeholder
  dimdat = 3
  numclust = 2
  TT = 3
  ## data = make_data(T)
  data = make_data(T, fac=0.5)
  ## par(mfrow=c(1,3))
  ## sapply(data, function(a)plot(a[,2:3], col=iris[,5], pch=16, cex=3))

  ## Initialize clusters
  la("~/repos/flowcy/flowcy")
  pie = init_pi(numclust, TT)
  pie = rbind(c(0.2,0.8), c(0.3, 0.7), c(0.7, 0.3))
  mu = init_mu(data, numclust, TT)
  sigma = init_sigma(data, numclust, TT, fac=1)
  n = nrow(data[[1]])
  s = 1E-6 #1E-4
  lam1 = 0 ## lam2 = 10

  ## Seeing what happens when drift EM
  la("~/repos/flowcy/flowcy")
  lam1 = 1000
  lam2 = 100
  results = driftem(data, mu, pie, niter=1000, sigma, T,
                    tol1 = 1E-10, tol2 = 1E-4, lam1, lam2, s, numclust)
  results_better
})

objects(results)
results$mulist[[19]]
myplot(results[[19]])
source("~/repos/flowcy/main/drift-helpers.R")
myplot(data, results$mulist[[results$final.iter]],
       sigma=results$sigmalist[[results$final.iter]], lam1=lam1, lam2=lam2,fac=20)
results$pielist

## Well, now it runs but something is seriously off. Find out what it is.
