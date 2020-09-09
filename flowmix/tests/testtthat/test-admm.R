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
    ##                     soft_threshC(a, b), times=10000)
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
    ##              expr     min       lq      mean  median      uq       max neval
    ##        rowSums(x) 151.614 167.3150 178.35363 168.771 179.096  3885.243 1e+05
    ##       rowSumsC(x)  26.741  41.8445  61.50820  55.570  60.090 52757.170 1e+05
    ##  rowSumsC_arma(x)  26.928  37.9450  60.56174  55.752  60.288  9104.738 1e+05



    ############################
    #### L2 ball projection ####
    ############################
    set.seed(0)
    mat = matrix(rnorm(300), nrow=100)
    C = 0.3
    pmat = projCmat(mat, C)
    pmatC = projCmatC(mat, C)
    expect_equal(pmat, pmatC)
    ## microbenchmark::microbenchmark(projCmat(mat, C),
    ##                                projCmatC(mat, C), times = 1000)
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


    ##############################################################################
    ## uw & Uz updates are simple matrix operations, so I'll leave them for now ##
    ##############################################################################
    ## uw <- uw_update(uw, rho, b1, wvec)
    ## Uw <- matrix(uw, nrow = p, byrow=FALSE)
    ## Uz <- Uz_update(Uz, rho, Xbeta1, Z)

    ## library(Rcpp)
    ## la('flowmix')
    ## sourceCpp("/home/shyun/repos/flowmix/flowmix/src/admm.cpp")
})


