context("Test M-step admm solver.")


test_that("Various substitutes give the same value", {


    #######################
    ## Soft thresholding ##
    #######################
    set.seed(0)
    a = rnorm(10)
    b = 1
    expect_equal(soft_thresh(a, b),
                 soft_threshC(a, b))

    ## microbenchmark::microbenchmark(
    ##                     soft_thresh(a, b),
    ##                     soft_threshC(a, b), times=100000)
    ## Unit: microseconds
    ##   expr  min    lq     mean median    uq    max neval
    ##   anew 5.17 5.862 6.484828  6.340 6.952 38.704 10000
    ##  anewC 1.93 2.164 2.498392  2.424 2.712 15.376 10000


    #################
    ## wvec update ##
    #################
    p = 10
    dimdat = 3
    set.seed(0)
    uw = rnorm(p * dimdat)
    b1 = rnorm(p * dimdat)
    rho = 1
    lambda = 1.5
    wvec <- wvec_update(b1, uw, lambda, rho)
    wvecC <- wvec_updateC(b1, uw, lambda, rho)
    expect_equal(wvec, wvecC)
    ## microbenchmark::microbenchmark(
    ##                     wvec_update(b1, uw, lambda, rho),
    ##                     wvec_updateC(b1, uw, lambda, rho), times = 10000)
    ## Unit: microseconds
    ##                               expr   min    lq     mean median     uq      max neval
    ##   wvec_update(b1, uw, lambda, rho) 6.761 7.528 8.826304 7.9475 9.1540 2582.323 10000
    ##  wvec_updateC(b1, uw, lambda, rho) 2.590 3.019 4.291510 3.2740 3.7135 8581.285 10000



    ##########################################
    ## Rowsums (probably already exists) #####
    ##########################################
    set.seed(1014)
    x <- matrix(sample(10000), 10)
    expect_equal(rowSums(x),
                 rowSumsC(x))
    expect_equal(rowSums(x),
                 rowSumsC_arma(x))
    microbenchmark::microbenchmark(rowSums(x),
                                   rowSumsC(x),
                                   rowSumsC_arma(x),
                                   times = 100000)
    ## microbenchmark::microbenchmark(rowSums(x),
    ##                                rowSumsC(x),
    ##                                rowSumsC_arma(x),
    ##                                times = 100000)

    ## Unit: microseconds
    ##              expr     min      lq      mean   median       uq        max neval
    ##        rowSums(x) 224.105 226.115 242.73947 227.0130 247.9290    955.312 1e+05
    ##       rowSumsC(x)  37.281  38.133  48.37031  38.9170  43.5825 280834.528 1e+05
    ##  rowSumsC_arma(x)  37.522  38.391  46.16909  39.1765  43.8940   7739.519 1e+05



    ############################
    #### L2 ball projection ####
    ############################
    set.seed(0)
    mat = matrix(rnorm(300), nrow=100)
    C = 0.3
    pmat = projCmat(mat, C)
    pmatC = projCmatC(mat, C)
    expect_equal(pmat, pmatC)
    sourceCpp("/home/shyun/repos/flowmix/flowmix/src/admm.cpp")
    microbenchmark::microbenchmark(projCmat(mat, C),
                                   projCmatC(mat, C), times = 10000)
    ## Unit: microseconds
    ##               expr    min     lq     mean  median     uq       max neval
    ##   projCmat(mat, C) 20.512 24.754 27.09897 26.2435 27.752   227.538 10000
    ##  projCmatC(mat, C)  6.465  9.394 12.59565 10.4140 11.207 17393.217 10000

    ## Unit: microseconds
    ##               expr    min      lq      mean  median     uq    max neval
    ##   projCmat(mat, C) 12.940 14.3415 16.468260 14.8975 15.659 65.842  1000
    ##  projCmatC(mat, C)  5.865  6.5805  7.702387  7.0240  7.549 41.808  1000


    ##############
    ## Z update ##
    ##############
    T = 100
    dimdat = 3
    set.seed(0)
    Xbeta1 = matrix(rnorm(TT * dimdat), nrow=TT)
    Uz = matrix(rnorm(TT * dimdat), nrow=TT)
    C = 0.1
    rho = 0.1
    Z <- Z_update(Xbeta1, Uz, C, rho, dimdat, TT)
    ZC <- Z_updateC(Xbeta1, Uz, C, rho, dimdat, TT)
    expect_equal(Z, ZC)
    ## microbenchmark::microbenchmark(
    ##                     Z_update(Xbeta1, Uz, C, rho, dimdat, TT),
    ##                     Z_updateC(Xbeta1, Uz, C, rho, dimdat, TT), times = 1000)
    ## Unit: microseconds
    ##                                       expr    min     lq     mean  median      uq    max neval
    ##   Z_update(Xbeta1, Uz, C, rho, dimdat, TT) 16.216 17.611 18.93144 18.2085 18.9280 53.860  1000
    ##  Z_updateC(Xbeta1, Uz, C, rho, dimdat, TT)  8.051  8.886  9.71677  9.4375 10.0495 26.498  1000

    ##############
    ## b update ##
    ##############

    set.seed(0)
    p = 30
    dimdat = 3
    TT = 100
    wvec = rnorm(p * dimdat)
    uw = rnorm(p * dimdat)
    rho = 0.1
    Z = rnorm(dimdat * TT)
    Uz = matrix(rnorm(dimdat * TT), ncol = dimdat)
    yvec = rnorm(dimdat * TT)
    D = matrix(rnorm(dimdat*(2*TT+p) * (p+1)*dimdat ), ncol = (p+1) * dimdat)
    DtDinv = solve(t(D) %*% D)
    N = 1000

    ## No speed improvement at all; in fact, slightly slower
    sourceCpp("/home/shyun/repos/flowmix/flowmix/src/admm.cpp")
    cvec3_el = as.numeric(t(Z - Uz/rho))
    bnewC = b_updateC(wvec, uw, rho, cvec3_el, yvec, DtDinv, D, N)##, dimdat, TT, p)
    bnew = b_update(wvec, uw, Z, Uz, rho, yvec, D, DtDinv, N)
    expect_equal(bnew, bnewC)

    ## microbenchmark::microbenchmark(
    ## {
    ##   cvec3_el = as.numeric(t(Z - Uz/rho));
    ##   bnewC = b_updateC(wvec, uw, rho, cvec3_el, yvec, DtDinv, D, N)
    ## },{
    ##   bnew = b_update(wvec, uw, Z, Uz, rho, yvec, D, DtDinv, N)
    ## }, times=100)
    ##Unit: microseconds
    ##     min       lq     mean   median       uq     max neval
    ## 178.126 198.6935 206.7051 200.7315 210.5055 381.198   100
    ## 151.098 165.5545 172.0651 171.3735 177.9080 212.127   100

    ##############################################################################
    ## uw & Uz updates are simple matrix operations, so I'll leave them for now ##
    ##############################################################################
    ## uw <- uw_update(uw, rho, b1, wvec)
    ## Uw <- matrix(uw, nrow = p, byrow=FALSE)
    ## Uz <- Uz_update(Uz, rho, Xbeta1, Z)
})


test_that("EM algorithm gives same results using ADMM solver or CVXR solver in M step.", {

  ## Generate data
  set.seed(0)
  TT = 300
  obj = generate_data_generic(nt = 3000, TT = TT, fac = 0.1, dimdat = 3)
  ylist = obj$ylist
  X = obj$X
  ## load("~/repos/cruisedat/export/MGL1704-hourly-only-binned.Rdata", verbose=TRUE)
  ## ylist = ybin_list
  ## countslist = biomass_list
  ## X = X %>% dplyr::select(-time, -lat, -lon) %>% as.matrix

  numclust = 10
  niter = 20
  mean_lambda=0.001
  pie_lambda=0.01
  print(mean_lambda)
  print(pie_lambda)
  library(profvis)
  ## ba()
  ## load_all()
  ## ## compileAttributes()
  ## load_all('/home/shyun/repos/flowmix/flowmix')

  niter = 1000
  la('flowmix')
  obj_rcpp = profvis({
    set.seed(2)
    res.admm = flowmix_once(ylist, X,
                            countslist=countslist,
                            numclust = numclust,
                            niter = niter,
                            mean_lambda = mean_lambda,
                            pie_lambda = pie_lambda,
                            verbose = TRUE,
                            maxdev = 0.5,
                            admm = TRUE,
                            admm_local_adapt_niter = 20, ## This spans rho=0.1 to
                            admm_rho = 10,
                            admm_err_rel = 1E-3,
                            admm_niter = 1E3,
                            rcpp = TRUE)
  })
  ## Method 2
  obj = profvis({
    set.seed(2)
    res.admm2 = flowmix_once(ylist, X,
                            countslist=countslist,
                             numclust = numclust,
                             niter = niter,
                             mean_lambda = mean_lambda,
                             pie_lambda = pie_lambda,
                             verbose = TRUE,
                             maxdev = 0.5,
                             admm = TRUE,
                            admm_local_adapt_niter = 20, ## This spans rho=0.1 to
                             admm_rho = 10,
                             admm_err_rel = 1E-3,
                             admm_niter = 1E3,
                             rcpp = FALSE)
  }, times = 1)

  ## Test
  expect_equal(res.admm, res.admm2)

  obj
  obj_rcpp
  save(obj, obj_rcpp, ylist, X, countslist, file=file.path("~/Desktop", "rcpp-comparison-profvis.Rdata"))
  load(file=file.path("~/Desktop", "rcpp-comparison-profvis.Rdata"), verbose = TRUE)
})
