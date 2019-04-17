##' Building a function to tune the number of GMM clusters.  The default for the
##' \code{Mclust()}  function  is to  use  BIC,  which  is  an estimate  of  the
##' marginalized    likelihood   value    in   the    same   data.     This   is
##' unsatisfactory. Instead, here,  we calculate the negative  log likelihood of
##' the data directly instead.
##' @param dat data table, of 3 columns
##' @param dat.new new data to evaluate likelihood at.
##' @param max.numclust Maximum number of clusters desired. Defaults to 80
##' @param  seed Seed  for fitting.   This is  because the  GMM fit  should use
##'   consistent initial values each time (for each G).
##' @param mc.cores number of multiple cores to use.
##' @param start.time Temporary
##' @return  A list object containing  (1) the negative log  likelihoods (2) the
##'   best number of clusters according to the
##' @export
##' @import mclust
tune <- function(dat, dat.new, max.numclust = 80, seed=NULL, mc.cores=1, start.time=NULL){

  ## If seed (for fitting) is not provided, set random seed
  if(is.null(seed)){ seed = as.numeric(Sys.time()) }
  Glist = 1:max.numclust

  ## Fit density of 1 through \code{max.numclust} clusters
  density.list = parallel::mclapply(Glist, function(G){
    printprogress(G, max.numclust, start.time=start.time) #Temporary! Get rid of this soon.
    set.seed(seed)
    dens <- densityMclust(dat, G=G, verbose=FALSE)
  }, mc.cores=mc.cores)

  ## Calculate negative log likelihood on \coe{dat.new}.
  pred.dens.list = sapply(1:max.numclust, function(G){
    dens = mclust::predict.densityMclust(density.list[[G]], dat.new[,])
    pred.dens = -sum(log(dens))
    return(pred.dens)
  })

  return(list(negloglik=pred.dens.list,
              best.numclust=which.min(pred.dens.list),
              dat=dat,
              dat.new=dat.new,
              max.numclust=max.numclust))
}


##' For  \code{newdata}  Find  the  best  mixture  proportions  to  the  mixture
##' components recovered in \code{obj}.
##' @param obj  An object of class \code{Mclust}, whose  mixture proportions are
##'   \code{G}.
##' @param newdata New data from which to obtain the mi
##' @param solver solver to use  for CVXR::solve(). Defaults to "SCS". For other
##'   options,                                                               see
##'   https://stackoverflow.com/questions/47586272/specify-solver-in-cvxr.
##' @return Vector of \code{G} fitted mixture (weights) proportions.
##' @import CVXR
##' @export
gmmreduce <- function(obj, newdata, solver="SCS"){

  all.dens = mclust::cdens(modelName = obj$modelName, data = newdata,
                           parameters = obj$parameters)

  ## Use CVXR to obtain best mixture weights
  p = obj$G
  prop <- CVXR::Variable(p)
  objective <- CVXR::Maximize(sum(log(all.dens %*% prop)))
  ## all.dens = all.dens[1:1000,]
  ## objective <- Maximize(sum(log(all.dens %*% prop)))
  constraint1 <- sum(prop)==1
  constraint2 <- prop >= 0
  problem <- CVXR::Problem(objective,
                           constraints = list(constraint1, constraint2))
  result <-solve(problem, solver=solver)
  sol <- as.numeric(result$getValue(prop))
  names(sol) = 1:p

  return(sol)
}


##' Produces function that takes in \code{newdata} (of the  same number of dimensions
##' as the original data, and returns the evaluated density at those points.
##' @param obj object of class |Mclust|
##' @param weights weights of each mixture component in \code{obj}.
##' @export
gmmreduce_density_fun <- function(obj, weights){

  stopifnot(obj$G == length(weights))

  ## This is a function that takes in new data, and calculates. Symbolically, it
  ## /is/ the resulting mixture density.
  return( function(dat){

    ## Basic checks
    stopifnot(ncol(obj$data) == ncol(dat))
    mns = (obj$parameters)$mean
    vars = (obj$parameters)$variance$sigma
    density.vals = lapply(1:obj$G, function(ii){
      a = mvtnorm::dmvnorm(dat, mns[,ii], vars[,,ii], log = FALSE)})
    weighted.densities = Map(function(x,y){x*y}, density.vals, weights)

    add <- function(x) Reduce("+", x)
    final.densities = add(Map(function(x,y){x*y}, density.vals, weights))

    ## Make colors between black and white
    colors = rescale_density(final.densities)

    return(list(each=weighted.densities,
                sum=final.densities,
                colors=colors))
  })
}


## Test case for the above
## set.seed(0)
## n = nrow(iris)
## isamp = sample(n, n/3)
## olddata = iris[isamp,1:3]
## newdata = iris[-isamp,1:3]
## all.dens = cdens(modelName = obj$modelName,
##                  data = newdata[c(1,30,40,43,47,50),],
##                  parameters = obj$parameters)
##  gmmreduce(obj, newdata)

## ##' Given newdata, \code{obj}
## sse_gmmreduce <- function(obj, sol){

## }



##' Process the  density so it  comes between  0 and 1,  and transforms it  to a
##' zero-one scale .
##' @param densities A vectors of density values (evaluated at
##'   desired data points).
##' @return Greyscale RGB colors.
##' @export
rescale_density <- function(densities){
  densities = densities/max(densities)
  f <- grDevices::colorRamp(c("white", "black"))
  colors <- grDevices::rgb(f(densities)/255)
  return(colors)
}

##' Obtain equally spaced grid in data range.  If \code{dat} has p columns, this
##' creates a numticks to the p power,  p-length grid locations.  Not to be used
##' if the data has more than 3 columns.
##' @param dat Data.
##' @param  numticks  The number  of  grid  values  to  take per  dimension  of
##'   dat. Defaults to 30.
##' @return A matrix of |numticks^p| rows and p columns, representing grid
##' points.
##' @export
get_grid <- function(dat, numticks=30){
  stopifnot(ncol(dat) <= 2)
  grids = apply(dat, 2, function(mycol){
    seq(from=min(mycol), to=max(mycol), length=numticks) })
  grids = expand.grid(data.frame(grids))
  return(grids)
}
