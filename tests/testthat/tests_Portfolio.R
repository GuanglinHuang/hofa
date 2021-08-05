library(testthat)

context("Portoflio.IC test")
test_that("Portoflio.IC test", {
  data("sp500")
  X = sp500[1:300,1:100]
  t = NROW(X)
  n = NCOL(X)
  result_ic = hofa::Portfolio.IC(X,r=20,Port_obj = "EU",Adjcov = "LI")

  expect_equal(tail(result_ic$obj,1),5.675809,tolerance = 10^-6)

})

context("Portoflio.PC test")
test_that("Portoflio.PC test", {
  data("sp500")
  X = sp500[1:300,1:100]
  t = NROW(X)
  n = NCOL(X)
  result_pc = hofa::Portfolio.PC(X,r=20,Port_obj = "EU",Adjcov = "LI")

  expect_equal(tail(result_pc$obj,1),5.865622,tolerance = 10^-6)

})
