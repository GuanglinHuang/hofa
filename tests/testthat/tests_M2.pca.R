library(testthat)

context("M2.pca test")
test_that("M2.pca test", {
  set.seed(1234)
  data("sp500")
  X = sp500[1:50,]
  t = NROW(X)
  n = NCOL(X)
  d = 4
  C = matrix(rnorm(n*d),n,d)
  result_pca = hofa::M2.pca(X, r = 2, method = "PCA")
  result_ppca = hofa::M2.pca(X,C,r = 2, method = "P-PCA")

  expect_equal(abs(result_pca$f[1,]),c(1.09396251,0.04828123),tolerance = 10^-6)
  expect_equal(abs(result_ppca$f[1,]),c(0.5313810,0.13233230),tolerance = 10^-6)

})
