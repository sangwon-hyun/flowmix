test_that("Test bootstrap and helpers", {

  ## Same-size checker
  dat1 = generate_data_generic(dimdat = 3)
  dat2 = generate_data_generic(dimdat = 3)
  expect_error(check_if_same_size(dat1$ylist, dat2$ylist), NA)

  ## Compatibility checker
  dat1 = generate_data_generic(dimdat = 2)
  res = flowmix_once(dat1$ylist, dat1$X, numclust=4, niter=100,
                     mean_lambda = 1E-3, prob_lambda = 1E-3)
  expect_error(check_compatible(ylist = dat1$ylist, res = res), NA)

})
