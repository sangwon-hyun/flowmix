context("Test eigendecompisition and supporting functions.")

## Test for lasso and cvx lasso
test_that("The conversion to eigenvalues behaves correctly.", {


  ## Create dummy data.
  set.seed(0)
  dat = matrix(rnorm(100),ncol=2)
  sigma = cov(dat)

  ## Test determinant.
  sigma_eig0 = eigendecomp_sigma_barebones(sigma)
  expect_equal(det_from_eig(sigma_eig0), det(sigma))

  ## Test inverse.
  expect_equal(sigma_inv_from_eig(sigma_eig0), solve(sigma))

  ## Test halving.
  sigma_half = sigma_half_from_eig(sigma_eig0)
  expect_equal(sigma_half %*% sigma_half, sigma)

  ## Test the multivariate normal density
  mu = apply(dat, 2, mean)
  sigma_eig = eigendecomp_sigma(sigma)
  sigma_eig0 = eigendecomp_sigma_barebones(sigma)
  expect_equal(dmvnorm_fast(y, mu, sigma_eig),
               mvtnorm::dmvnorm(y, mu, sigma))


  ## Test that the eigenspeed results are the same as before
  set.seed(1)
  dat = generate_data_generic(TT=500, p=5, nt=2000)
  ylist = dat$ylist
  X = dat$X
  set.seed(0)
  obj1 = covarem(ylist, X, numclust=4, eigenspeed=TRUE, niter=5, nrep=1, verbose=TRUE)
  set.seed(0)
  obj2 = covarem(ylist, X, numclust=4, eigenspeed=FALSE, niter=5, nrep=1, verbose=TRUE,
                 faster_mvn=TRUE)
  expect_true(max(max(obj1$alpha-obj2$alpha)) < 1E-4) ## I don't get an exact match.

})
