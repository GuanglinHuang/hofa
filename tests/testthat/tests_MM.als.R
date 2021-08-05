library(testthat)

context("M3.als test")
test_that("M3.als test", {
  data("sp500")
  X = sp500[1:300,1:100]
  t = NROW(X)
  n = NCOL(X)
  result_als = hofa::M3.als(X,gamma = c(0,1),rh = 1,rg = 1)

  expect_equal(abs(result_als$f[1,]),c(1.1592216, 0.1614519),tolerance = 10^-6)

})


context("M4.als test")
test_that("M4.als test", {
  data("sp500")
  X = sp500[1:300,1:100]
  t = NROW(X)
  n = NCOL(X)
  result_als = hofa::M4.als(X,gamma = c(0,0,1),rh = 1,rg = 1)

  expect_equal(abs(result_als$f[1,]),c(1.12940367, 0.06014121),tolerance = 10^-6)

})
