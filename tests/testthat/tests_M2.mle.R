library(testthat)
set.seed(1)

test_that("M2.mle test", {
  data("sp500")
  X = sp500[1:300,1:100]
  t = NROW(X)
  n = NCOL(X)

  result_ml = hofa::M2.mle(X,r =2,method = "ML")
  result_qml = hofa::M2.mle(X,r =2,method = "QML")
  result_gls = hofa::M2.mle(X,r =2,method = "ML-GLS")
  result_ite = hofa::M2.mle(X,r =2,method = "ML-ITE")
  result_em = hofa::M2.mle(X,r =2,method = "ML-EM")

  abs(result_ml$f[1,])
  abs(result_qml$f[1,])
  abs(result_gls$f[1,])
  abs(result_ite$f[1,])
  abs(result_em$f[1,])

  expect_equal(abs(result_ml$f[1,]),c(0.8341876 , 0.7186874),tolerance = 10^-6)
  expect_equal(abs(result_qml$f[1,]),c(0.8341876 , 0.7186874),tolerance = 10^-6)
  expect_equal(abs(result_gls$f[1,]),c(0.8328814 , 0.6959627),tolerance = 10^-6)
  expect_equal(abs(result_ite$f[1,]),c(0.8286102 , 0.7043564),tolerance = 10^-6)
  expect_equal(abs(result_em$f[1,]),c(0.30867413, 0.04451294),tolerance = 10^-6)

})
