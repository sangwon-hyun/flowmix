test_that("New m step ADMM solver", {



  sourceCpp("~/repos/flowmix/flowmix/src/estep.cpp")
  m = matrix(rnorm(90), ncol=3);
  ss(m, c(1, 1, 5, 12, 22))

  ## Generate / Load data
  load("~/repos/cruisedat/export/MGL1704-hourly-only-binned.Rdata", verbose=TRUE)
  ylist = ybin_list
  countslist = biomass_list
  X = X %>% dplyr::select(-time, -lat, -lon) %>% as.matrix
  numclust = 10
  niter = 100
  maxdev = 0.5
  mean_lam = 0.001
  prob_lam = 0.01
  niter = 20
  admm_rho = 0.01
  la('flowmix')
  speed.new.aggtype4 = profvis::profvis({
    set.seed(2)
    obj.new.aggtype3 = flowmix_once(ylist=ylist,
                           countslist = countslist,
                           X=X, numclust=numclust,
                           prob_lambda = prob_lam,
                           mean_lambda = mean_lam,
                           maxdev = maxdev,##0.5,
                           niter = niter,
                           verbose = TRUE,
                           admm_rho = admm_rho)
  })
  speed.new.aggtype1
  ## speed.new.aggtype2
  speed.new.aggtype4
  range(obj.new.aggtype1$prob -   obj.new.aggtype2$prob)
  range(obj.new.aggtype1$prob -   obj.new.aggtype3$prob)
  speed.new2
  ## save(obj.new, speed.new, file="~/Desktop/slow-obj.Rdata")
  ## save(obj.new, speed.new, file="~/Desktop/fast-obj.Rdata")
  ## load(file="~/Desktop/slow-obj.Rdata")
  ## load(file="~/Desktop/fast-obj.Rdata")

  speed.old = profvis::profvis({
    set.seed(2)
    obj.old = flowmix_once(ylist=ylist,
                           countslist = countslist,
                           X=X, numclust=numclust,
                           prob_lambda = prob_lam,
                           mean_lambda=mean_lam,
                           new_beta_mstep = FALSE,
                           maxdev = maxdev, ##0.5
                           niter = niter,
                           verbose=TRUE,
                           admm_rho = admm_rho,
                           rcpp = TRUE
                           ## ## The previous settings, exactly!
                           ## admm_local_adapt_niter = 20,
                           ## admm_rho = 10,
                           ## admm_err_rel = 1E-3,
                           ## admm_niter = 1E3
                           )
  })

  speed.new
  speed.old
  obj.new$beta
  save(speed.old, obj.old, speed.new, obj.new,
       file=file.path("~/Desktop/admm-speed-rho10.Rdata"))

  speed.old
  speed.new
  obj.new$niter
  obj.old$niter

  obj.new %>% objects()
  plot(obj.new$prob, obj.old$prob)
  plot(obj.new$mn,obj.old$mn, pch=16)##, cex=.2)
  abline(0,1)
  par(mfrow=c(1,2))
  matplot(obj.new$mn[,1,, drop=TRUE])
  matplot(obj.old$mn[,1,, drop=TRUE])
  abline(0, 1)
  par(mfrow=c(2,5))
  for(iclust in 1:10){
  plot(obj.old$beta[[iclust]], obj.new$beta[[iclust]])
  abline(0,1)
  }

  dim(obj.new$mn)
  obj.new %>% objects()
  plot(obj.new$objectives)
  lines(obj.old$objectives)






  ## Load things right before mstep is to be called
  load(file = file.path("~/Desktop", "ADMM-test.Rdata"), verbose=TRUE)

  ## Source in all helpers
  Sys.setenv("PKG_LIBS" = "-llapack")
  Rcpp::compileAttributes("~/repos/flowmix/flowmix")
  library(Rcpp)
  la('flowmix')
  sourceCpp("~/repos/flowmix/flowmix/src/syl.cpp")
  maxdev = 1000000
  ## lambda = 0

  ## beta = beta_update_new(schurA, schurB, syl_C)
  ## Rcpp::compileAttributes("~/repos/flowmix/flowmix")
  ## obj.new = profvis::profvis({
  ##   ## time.after = microbenchmark::microbenchmark({
      res.beta = Mstep_beta_admm_new(resp, ylist, X,
                                     mean_lambda = lambda,
                                     first_iter = TRUE,
                                     sigma = sigma, maxdev = maxdev,
                                     rho = 0.1)
    ## }, times=10)
  })

 ## obj.old = profvis::profvis({
   time.before = microbenchmark::microbenchmark({
    res.beta.old = Mstep_beta_admm(resp, ylist, X,
                                 mean_lambda = lambda,
                                 first_iter = TRUE,
                                 sigma = sigma, maxdev = maxdev,
                                 rho = 0.1,
                                 rcpp = FALSE)
    }, times=10)
  ## })

  obj.new
  obj.old

  plot(abs(res.beta$betas[[1]]),
       abs(res.beta.old$betas[[1]]),
       log = "xy")
  abline(0,1)

  obj.old.rcpp = profvis::profvis({
    res.beta.old = Mstep_beta_admm(resp, ylist, X,
                                 mean_lambda = lambda,
                                 first_iter = TRUE,
                                 sigma = sigma, maxdev = maxdev,
                                 rho = 0.1,
                                 rcpp = TRUE)
  })

  obj.new
  obj.old
  obj.old.rcpp

  (res.beta.old$beta[[1]] - res.beta$beta[[1]]) %>% range()

  plot(y=rbind(beta0, beta$beta) %>% abs(), ylab = "new",
       x=res.beta$beta[[iclust]] %>% abs(), xlab = "old", log="xy")
  abline(0,1)

  fit_new = objective_per_cluster_temp(rbind(beta0, beta$beta),
                                       ylist, Xa, resp,
                                       lambda, N, dimdat,
                                       iclust, sigma, iter, zerothresh)
  fit_old = objective_per_cluster_temp(res.beta$beta[[iclust]],
                                       ylist, Xa, resp,
                                       lambda, N, dimdat,
                                       iclust, sigma, iter, zerothresh)
  fit_new
  fit_old

  ## Also check that a single beta update is identical

  ## New
  ## Initialize some variables
  load(file = file.path("~/Desktop", "ADMM-test.Rdata"), verbose=TRUE)
  Z = matrix(0, nrow = TT, ncol = dimdat)
  W = matrix(0, nrow = p, ncol = dimdat)
  U = matrix(0, nrow = TT + p, ncol = dimdat)
  beta = matrix(0, nrow = p, ncol = dimdat)
  beta0 = rep(0, dimdat)
  zrows = 1:TT
  wrows = TT + (1:p)
  Xa = cbind(1, X)
  zerothresh = 1E-6
  iclust = 2

  ## Some intermediate quantities
  resp.sum.thisclust = sum(resp.sum[,iclust])
  ybar = Reduce("+", Map(function(y, myresp){
    ybar = (myresp[,iclust, drop = TRUE] ) %*% y
  }, ylist, resp)) / resp.sum.thisclust
  ycentered_list = Map(function(y, myresp){
    return(colSums(myresp[,iclust] * sweep(y, 2, ybar)))
  }, ylist, resp)
  ycentered = do.call(cbind, ycentered_list)
  Xtilde = colSums(resp.sum[,iclust] * X) / resp.sum.thisclust
  Xcentered = sweep(X, 2, Xtilde)
  Xcentered %>% dim()
  D = diag(resp.sum[,iclust])

  source("~/repos/flowmix/flowmix/R/mstep_new.R")
  rho = 0.1
  beta_new = beta_update(beta, ylist, rho, sigma, X, Xinv, Xaug, iclust, U, Z,
                         W, ycentered, Xcentered, D)
  beta0_new = intercept(resp, ylist, beta_new, X, N, iclust)
  beta_new = rbind(beta0_new,
                   beta_new)

 ## Old

  ####################
  load(file = file.path("~/Desktop", "ADMM-test.Rdata"), verbose=TRUE)
  TT = length(ylist)
  p = ncol(X)
  numclust = ncol(resp[[1]])
  dimdat = ncol(ylist[[1]])
  ntlist = sapply(ylist, nrow)
  ## N = sum(ntlist) ## OLD BUG!!
  resp.sum = t(sapply(resp, colSums)) ## (T x numclust)
  resp.sum = as.matrix(resp.sum)
  N = sum(resp.sum) ## NEW (make more efficient, later)

  Xa = cbind(1, X)
  intercept_inds = ((1:dimdat) - 1)*(p+1) + 1
  X0 = lapply(1:TT, function(tt){ diag(rep(1,dimdat)) %x% t(c(0, X[tt,,drop=TRUE]))})
  X0 = do.call(rbind, X0)
  tX = t(X)
  I_aug = make_I_aug(p, dimdat, intercept_inds)

  ## Form tilde objects for b update. Only do once!
  manip_obj = manip(ylist, Xa, resp, sigma, numclust,
                    sigma_eig_by_clust = sigma_eig_by_clust,
                    first_iter = TRUE)
  Xtildes = manip_obj$Xtildes
  yvecs = manip_obj$yvecs
  Xtilde = Xtildes[[iclust]]
  yvec = yvecs[[iclust]]

  Dfirst = sqrt(1/(2*N)) * Xtilde
  Drest = rbind(sqrt(rho/2) * I_aug,
                sqrt(rho/2) * X0)
  D = rbind(Dfirst, Drest)
  DtD = crossprod(Dfirst, Dfirst) + crossprod(Drest, Drest)
  DtDinv = chol2inv(chol(DtD))
  Dobj = DtDinv %*% t(D)

  ###############################
  ## Initialize the variables ###
  ###############################
  C = maxdev
  converge = FALSE
  Z = matrix(0, nrow = TT, ncol = dimdat)
  wvec = rep(0, p * dimdat)
  uw  = rep(0, p * dimdat)
  Uz = matrix(0, nrow = TT, ncol = dimdat)

  b = b_update(wvec, uw, Z, Uz, rho, yvec, D, DtDinv, N)
  b1 = b[-intercept_inds]
  b0 = b[intercept_inds]
  beta1 = matrix(b1, nrow = p)
  beta_old = rbind(b0, beta1)

  par(mfrow=c(2,4))
  logg = "xy"
  plot(abs(beta_old), abs(beta_new), log=logg)
  abline(0,1)
  plot(abs(beta_old)[,1], abs(beta_new)[,1], log=logg)
  abline(0,1)
  plot(abs(beta_old)[,2], abs(beta_new)[,2], log=logg)
  abline(0,1)
  plot(abs(beta_old)[,3], abs(beta_new)[,3], log=logg)
  abline(0,1)
  plot(abs(cbind(X)%*% beta_old[-1,]), abs(cbind(X) %*% beta_new[-1,]), log=logg)
  abline(0,1)
  for(idim in 1:3){
    plot(abs(cbind(X)%*% beta_old[-1,])[,idim], abs(cbind(X) %*% beta_new[-1,])[,idim], log=logg,
         main=idim)
    abline(0,1)
  }

 ## Calculate the objective/
  fit_new = objective_per_cluster_temp(beta_new,
                                       ylist, Xa, resp,
                                       lambda, N, dimdat,
                                       iclust, sigma, iter, zerothresh)
  fit_old = objective_per_cluster_temp(beta_old,
                                       ylist, Xa, resp,
                                       lambda, N, dimdat,
                                       iclust, sigma, iter, zerothresh)

  load(file="~/Desktop/new-mstep-beta.Rdata")
  load(file="~/Desktop/beta-mstep-former.Rdata")


  beta0_old = beta0
  beta_old = beta

  plot(x=rbind(beta0, beta), y=beta_old)##rbind(beta0_old, beta_old))
  abline(0,1)
  length(beta)
  length(beta_old)

  plot(res.beta$betas[[iclust]], rbind(beta0, beta))
  abline(0,1)
  diffs = res.beta$betas[[iclust]] - rbind(beta0, beta)
  inds = which(abs(diffs) < 1E-2)##, arr.ind=TRUE)
  cols = rep("black", length(diffs))
  cols[inds] = "red"
  plot(res.beta$betas[[iclust]], rbind(beta0, beta), col=cols)
  abline(0,1)
  inds %>% length()
  diffs %>% length()
  round(diffs,2)

  ## This is weird because diam_mid is all correct, and all


})
