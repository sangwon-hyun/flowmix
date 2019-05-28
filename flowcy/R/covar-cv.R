##' Make all CV indices.
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
  ## objective(ylist.test, res.train)

  ## Assign mn and pie
  pred = predict.covarem(res.train)

  ## Calculate objective (penalized likelihood)
  objective_overall_cov(aperm(pred$newmn, c(1,3,2)),
                        pred$newpie, pred$sigma, ylist.test,
                        pie_lambda=pie_lambda,
                        mean_lambda=mean_lambda,
                        alpha=res.train$alphalist[[res.train$final.iter]],
                        beta=res.train$betalist[[res.train$final.iter]])
}

##' CV wrapper for covarem().
##' @param lambdas_means Regularilization parameter for mean.
##' @param lambdas_pies Regularizization paraemter for pie.
##' @param nsplit Number of CV splits. Defaults to 5.
##' @param ... default arguments to covarem().
cv.covarem <- function(ylist=ylist, X=X, mean_lambdas, pie_lambdas, nsplit=5, ...){

  mysplits = cvsplit(ylist, nsplit=5)
  meanmat = matrix(NA, ncol=length(mean_lambdas), nrow=length(pie_lambdas))
  sdmat = matrix(NA, ncol=length(mean_lambdas), nrow=length(pie_lambdas))

  ## Main CV double loop (nothing special)
  for(imean in 1:length(mean_lambdas)){
    for(ipie in 1:length(pie_lambdas)){

      ## Assign mean and pie
      mean_lambda = mean_lambdas[imean]
      pie_lambda = pie_lambdas[ipie]
      print(c(mean_lambda, pie_lambda))

      ## Do CV.
      res = get_cv_score(ylist, X, mysplits, nsplit,
                   ## Other arguments to covarem go here.
                   mean_lambda=mean_lambda,
                   pie_lambda=pie_lambda,
                   ...)


      meanmat[imean, ipie] = res$mean
      sdmat[imean, ipie] = sd(res$all)
    }
  }
  return(list(meanmat=meanmat, sdmat=sdmat))
}
