## Generate data
datobj = generate_data_generic(p=5, TT=30, fac=0.05, nt=100)
set.seed(0)
ylist = datobj$ylist
X = datobj$X
numclust = 4

########################################################################
## Run algorithm one time, for a fixed pair of regularization
## parameters. Change |nrep| for multiple restarts.
########################################################################
set.seed(0)
res = covarem(ylist, X, numclust = numclust, niter = 500,
              mean_lambda = 0.01,
              pie_lambda = 10,
              verbose = TRUE,
              maxdev = NULL,
              nrep = 1)
print(res)
print(res$beta)
print(res$alpha)
plot_obj = fancyplot(res)
plot_obj
