## Synopsis: contains the main CV wrapper, and all helper functions related to
## cross-validations.

##' CV wrapper for covarem().
##' @param nsplit Number of CV splits. Defaults to 5.
##' @param ... default arguments to covarem().
##' @return List containing (1) the set of coefficients
cv.covarem <- function(ylist=ylist, X=X, mean_lambdas, pie_lambdas,
                       nsplit=5,
                      ...){
  ## Basic checks
  stopifnot(length(mean_lambdas) == length(pie_lambdas))
  gridsize = length(mean_lambdas)

  ## Create CV split indices
  mysplits = cvsplit(ylist, nsplit=5)

  ## Empty containers
  meanmat = sdmat = matrix(NA, ncol=gridsize, nrow=gridsize)
  cvscoremat = array(NA, dim=c(gridsize,gridsize,5))

  reslist = list()
  for(ii in gridsize:1){
    for(jj in gridsize:1){
      print(c(ii,jj))

      ## Set warm starts.
      if(jj == gridsize){
        if(ii == gridsize)  warmstart_mn = NULL
        if(ii != gridsize)  warmstart_mn = reslist[[paste0(ii+1, "-", jj)]]$mn
      } else {
        warmstart_mn = reslist[[paste0(ii, "-", jj+1)]]$mn
      }

      ## Do CV to get score
      cvres = get_cv_score(ylist, X, mysplits, nsplit,
                           mean_lambda=mean_lambdas[ii],
                           pie_lambda=pie_lambdas[jj],
                           mn=warmstart_mn,
                           ...)

      ## Get the fitted results on the entire data
      res = covarem(ylist=ylist, X=X,
                    mean_lambda=mean_lambdas[ii],
                    pie_lambda=pie_lambdas[jj],
                    mn=warmstart_mn, ...)

      ## Store the results
      reslist[[paste0(ii, "-", jj)]] = res

      ## Record CV score
      meanmat[ii, jj] = cvres$mean
      sdmat[ii, jj] = sd(cvres$all)
      cvscoremat[ii,jj,]=cvres$all
    }
  }
  return(list(reslist=reslist,
              meanmat=meanmat,
              sdmat=sdmat,
              cvscoremat=cvscoremat))
}


##' Makes all CV indices.
##' @param ylist ylist.
##' @param nsplit Defaults to 5.
##' @param return TT-length list of indices to use for CV.
cvsplit<- function(ylist, nsplit=5){
  TT = length(ylist)
  all.cv.inds = lapply(1:TT, function(tt){
    flds <- caret::createFolds(ylist[[tt]][,1], k = nsplit, list = TRUE, returnTrain = FALSE)
    names(flds)[1] <- "train"
    return(flds)
  })
}

##' Getting cross-validated test likelihood, averaged over train/test splits.
##' @param splits TT-lengthed list of indices.
##' @param ylist List of responses.
##' @param ... arguments to covarem
##' @return Cross validated test likelihood, scalar valued.
get_cv_score <- function(ylist, X, splits, nsplit,...){
  ## stopifnot(length(splits[[1]])!=nsplit) ## good check but only works if TT>1

  ## Cycle through splits, and calculate CV scroe
  all.scores = sapply(1:nsplit, function(test.isplit){
    get_cv_score_onesplit(test.isplit, splits, ylist, X,...)
  })
  return(list(mean=mean(all.scores), all=all.scores))
}

##' Inner  function  for calculating  cross-validated  test  likelihood for  one
##' train/test split.   Specifically, it takes  the set of splitted  indices for
##' all times in 1:TT, trains on (1:nsplit)[-isplit], tests on splits
##' @param splits TT-lengthed list of indices.
##' @param test.split The split (out of 1:nsplit) to use for test.
##' @param ylist List of responses.
##' @param ... arguments to covarem
##' @return One split's test likelihood.
get_cv_score_onesplit <- function(test.isplit, splits, ylist, X, pie_lambda, mean_lambda,...){

  TT = length(ylist)

  ## Obtain train and test data according to test split
  ylist.test = lapply(1:TT, function(tt){
    ind = splits[[tt]][[test.isplit]]
    ylist[[tt]][ind,]
  })
  ylist.train = lapply(1:TT, function(tt){
    ind = splits[[tt]][[test.isplit]]
    ylist[[tt]][-ind,]
  })

  ## Run algorithm on train data, evaluate test data.
  res.train = covarem(ylist.train, X, pie_lambda=pie_lambda, mean_lambda=mean_lambda,
                      ...)

  ## Assign mn and pie
  pred = predict.covarem(res.train)
  stopifnot(all(pred$newpie)>=0)

  ## Calculate objective (penalized likelihood)
  objective_overall_cov(aperm(pred$newmn, c(1,3,2)),
                        pred$newpie, pred$sigma, ylist.test,
                        pie_lambda=pie_lambda,
                        mean_lambda=mean_lambda,
                        alpha=res.train$alphalist[[res.train$final.iter]],
                        beta=res.train$betalist[[res.train$final.iter]])
}


##' Create a list containing the candidate regularization parameters for EM.
make_lambdas <- function(ylist, X, numclust, cv.grid.size=5){
  res0 = get_param0(ylist, numclust)
  max_lambda_beta = lambda_beta_max(res0$alpha0, res0$beta0, res0$sigmalist0,
                                    numclust, X, res0$resplist0, ylist)
  lambdas_beta = c(0,exp(seq(from=0, to=log(max_lambda_beta),
                             length=cv.grid.size)))
  max_lambda_alpha = lambda_alpha_max(res0$alpha0, res0$beta0, res0$sigmalist0,
                                      numclust, X, ylist)
  lambdas_alpha = c(0,exp(seq(from=0, to=log(max_lambda_alpha), length=cv.grid.size)))
  return(list(lambdas_beta=lambdas_beta,
              lambdas_alpha=lambdas_alpha))
}
