####################
## Generate data ###
####################
datobj = generate_data_generic(p=5, TT=300, fac=1, nt=2000)
set.seed(0)
ylist = datobj$ylist
X = datobj$X
numclust = 4

########################################################################
## Run algorithm one time, for a fixed pair of regularization
## parameters. Change |nrep| for multiple restarts.
########################################################################
set.seed(0)
res = flowmix_once(ylist, X, numclust = numclust, niter = 5,
              mean_lambda = 0.01,
              prob_lambda = 10,
              verbose = TRUE,
              maxdev = 0.5)
print(res)

########################################################################
## 2. Cross validate over a 2d grid of regularization parameters.
########################################################################
gridsize = 12
destin = "~" ## Where to save the CV results to

## Get range of lambdas
maxres = get_max_lambda(destin = destin,
                        ylist = ylist,
                        countslist = NULL,
                        X = X, numclust = numclust, verbose = TRUE,
                        max_lambda_alpha = 2, max_lambda_beta = 2, maxdev = 0.5)
prob_lambdas = seq(from  =  0, to = maxres$alpha, length = gridsize)
mean_lambdas = seq(from = 0, to = maxres$beta, length = gridsize)


## The rest of the code is under construction; the cv function below needs to be
## modified to run on a single computer.

## ## Run the parallel CV. Define cl, multicore.cv, etc.
## cv.flowmix(ylist = ylist, X = X, numclust = numclust,
##                     prob_lambdas = prob_lambdas,
##                     mean_lambdas = mean_lambdas,
##                     maxdev = 0.5,
##                     ## Options for parallelizing
##                     gridsize = gridsize,
##                     nsplit = 5,
##                     numfork = 12,
##                     verbose = TRUE,
##                     nrep = 10,
##                     destin = "~",
##                     warm_start = FALSE,
##                     cl = NULL)

## aggregateres_df(gridsize, destin)
## get_optimal_info(isim, outputdir = destin,
##                  gridsize = 12, excludecorner = FALSE)

## ##################################################
## ## Run the algorithm with binning *on the fly* ###
## ##################################################
## ylist_collapsed = do.call(rbind, ylist)
## dimdat = ncol(ylist[[1]])
## ranges = lapply(1:dimdat, function(ii) range(ylist_collapsed[,ii])) ## Get overall range.
## dat.grid = make_grid(ranges, gridsize = dat.gridsize)
## flowmix_once(ylist, X, numclust=4, verbose = TRUE,
##              prob_lambda = .1,
##              mean_lambda = .1,
##              manual.bin = TRUE,
##              manual.grid = dat.grid)
