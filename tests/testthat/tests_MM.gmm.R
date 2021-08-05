library(testthat)

context("M2.gmm test")
test_that("M2.gmm test", {
  data("sp500")
  X = sp500[1:300,1:100]
  t = NROW(X)
  n = NCOL(X)
  result_gmm = hofa::M2.gmm(X,r = 2,kappa = 0.0,initial = "MLE")

  expect_equal(abs(result_gmm$f[1,]),c(1.09585954, 0.04196591),tolerance = 10^-6)

})

context("M3.gmm test")
test_that("M3.gmm test", {
  data("sp500")
  X = sp500[1:300,1:100]
  t = NROW(X)
  n = NCOL(X)
  result_gmm = hofa::M3.gmm(X,r = 2,kappa = 0.0,initial = "MLE")

  expect_equal(abs(result_gmm$f[1,]),c(1.09481362, 0.05183387),tolerance = 10^-6)

})

