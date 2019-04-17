context("Test optimization for the drifting clusters problem")

test_that("Initialization works well.", {
  data(iris)
  dimdat = 3
  datapoints = iris[,1:dimdat]
  numclust = 2
  T = 3
  pie = init_pi(numclust, T)
  mu = init_mu(data, numclust, T)
  sigma = init_sigma(data, numclust, T)
  testthat::expect_equal(dim(pie), c(T, numclust))
  testthat::expect_equal(dim(mu), c(T, numclust, dimdat))
  testthat::expect_equal(dim(sigma), c(T, numclust, dimdat, dimdat))
})
