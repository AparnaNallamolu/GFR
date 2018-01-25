library(Prueba)
library(BGLR)
context("Prueba BGLR")

load("HTP_Irrigated.RData", verbose = T)
ETA <- list(X<-)
test_that("BFDA works like BGLR:", {

  expect_error(1 / "a", "To realice Fold Cross-validation BGFRA requieres y param like a data.frame with Line and Env specified on it")
  expect_equal(BFDA(), BGLR())
  expect_equal(BFDA(), BGLR())
})
#
# test_that("str_length of factor is length of level", {
#   expect_equal(str_length(factor("a")), 1)
#   expect_equal(str_length(factor("ab")), 2)
#   expect_equal(str_length(factor("abc")), 3)
# })
#
# test_that("str_length of missing is missing", {
#   expect_equal(str_length(NA), NA_integer_)
#   expect_equal(str_length(c(NA, 1)), c(NA, 1))
#   expect_equal(str_length("NA"), 2)
# })
