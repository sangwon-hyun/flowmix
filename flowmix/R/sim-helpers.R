########################################################################################
## All just convenience functions, temporarily. Nothing in this file will be retained ##
########################################################################################

## Let me reproduce the data here:
make_data75 <- function(){
  ## Load the hourly level X
  datadir = "~/repos/cruisedat/export"
  load(file.path(datadir, "MGL1704-hourly-only-binned.Rdata"))

  ## For back-compatibility; X now doesn't have any NA par values.
  if(any(is.na(X[,"par"]))){
    exclude_ind = which(is.na(X[,"par"]))
    X = X[-exclude_ind,]
  }

  ## Create the TF bases.
  TT = nrow(X)
  ind = rep(0, TT)
  north = which(X[,"lat"] > 37)
  regions = sapply(1:3, function(ii){
    endpt = c(0,range(north),TT)[ii:(ii+1)]
    (endpt[1]+1):endpt[2]
  })

  load(file.path(datadir, "MGL1704-hourly-only-binned-1d-diam.Rdata"))
  print("Hourly time resolution with bases and lagged PAR")
  print("Only sss, sst, par are being used.")
  X = X[,which(colnames(X) %in% c("sss", "sst", "par", "sss_cruise", "sst_cruise"))]
  print("Artificial covariates are being added.")

   ## Create the lagged sunlight variable
  lags = c(0,3,6,9,12)
  par = scale(X[,"par"])
  par = par - min(par)
  par = par / max(par)
  parlist = lapply(lags, function(lag)lagpad(par, lag))

  ## Make the additional columns
  dat = do.call(cbind, parlist)
  colnames(dat) = paste0("p", 0:4)

  ## Combine them with X
  X = X[,-which(colnames(X) == "par")]
  X = cbind(X, dat)
  print(colnames(X))
  X = scale(X)

  ## Add bases
  bases0 = lapply(regions, function(reg){
    vec = rep(0, TT)
    vec[reg] = 1
    vec
  })
  bases = do.call(cbind, bases0)
  ## bases = do.call(cbind, bases0[2:3]) ## Omitting one basis.
  ## colnames(bases) = c("b2", "b3")
  colnames(bases) = c("b1", "b2", "b3")
  X = cbind(X, bases)

  ## Rid of NA rows.
  na.rows = which(apply(dat, 1, function(myrow)any(is.na(myrow))))
  if(length(na.rows) > 0){
    X = X[-na.rows,]
    ylist = ylist[-na.rows]
    countslist = countslist[-na.rows]
    ylist = lapply(ylist, cbind)
  }

  maxdev = 0.1155245
  return(list(ylist=ylist,
              X=X,
              countslist=countslist))
}





##' Adding noise to covariates.
##'
##' @param Xorig Original covariate matrix
##' @param sigmalist Numeric vector of Gaussian noise standard deviations.
##' @param noise_ii Index of \code{sigmalist} to use.
##'
##' @return A noisy version of \code{Xorig}.
##'
##' @export
add_noise <- function(Xorig, sigmalist, noise_ii){

  X = Xorig

  if(noise_ii > 0){
    ## Add noise to the sunlight covariate
    X[,"par"] = Xorig[,"par"] + rnorm(TT, 0, sigmalist[noise_ii])

    ## Add noise to the spurious covariates by the same amount as well.
    for(icol in 3:ncol(X)){
      X[,icol] = Xorig[,icol] + rnorm(TT, 0, sigmalist[noise_ii])
    }
  }
  return(X)
}


##' Creates a directory \code{destin}, if it doesn't already exist.
##' @export
create_destin <- function(destin){
  if(!dir.exists(destin)){
    dir.create(destin, recursive = TRUE)
    cat("Creating destin: ", destin, fill=TRUE)
  } else {
    cat("All output goes out to destin: ", destin, fill = TRUE)
  }
}




##' Helper to divide up the jobs in \code{iimat} into a total of
##' \code{arraynum_max} jobs. The purpose is to divide the jobs, in order to run
##' this on a server.
##'
##' @param arraynum_max Maximum SLURM array number.
##' @param iimat matrix whose rows contain CV job indices.
##'
##' @export
make_iilist <- function(arraynum_max, iimat){
  iimax = nrow(iimat)
  if(arraynum_max > iimax){
    iilist = lapply(1:iimax, function(a)a)
  } else {
    ends = round(seq(from=0, to=iimax, length=arraynum_max+1))
    iilist = Map(function(a,b){ (a+1):b}, ends[-length(ends)], ends[-1])
    stopifnot(length(unlist(iilist)) == nrow(iimat))
  }
  stopifnot(length(unlist(iilist)) == nrow(iimat))
  stopifnot(all(sort(unique(unlist(iilist))) == sort(unlist(iilist))))
  return(iilist)
}
