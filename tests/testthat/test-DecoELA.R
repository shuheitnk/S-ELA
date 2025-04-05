test_that("DecoELA correctly calculates features", {
  X <- matrix(runif(200), nrow = 100)
  Y <- matrix(runif(200), nrow = 100)
  deco_features <- DecoELA(X = X, Y = Y, aggregate = TRUE, H = 5, scalar_func = "weightedsum", set_name = "ela_meta")
  expect_is(deco_features, "list")
  expect_true(all(grepl("^deco\\.", names(deco_features))))
})

test_that("DecoELA correctly calculates features (ela_distr)", {
  X <- matrix(runif(200), nrow = 100)
  Y <- matrix(runif(200), nrow = 100)
  deco_features <- DecoELA(X = X, Y = Y, aggregate = TRUE, H = 5, scalar_func = "weightedsum", set_name = "ela_distr")
  expect_is(deco_features, "list")
  expect_true(all(grepl("^deco\\.", names(deco_features))))
})

test_that("DecoELA correctly calculates features (ela_meta, basic)", {
  X <- matrix(runif(200), nrow = 100)
  Y <- matrix(runif(200), nrow = 100)
  deco_features <- DecoELA(X = X, Y = Y, aggregate = TRUE, H = 5, scalar_func = "weightedsum", set_name = "basic")
  expect_is(deco_features, "list")
  expect_true(all(grepl("^deco\\.", names(deco_features))))
})

test_that("DecoELA correctly calculates features (disp)", {
  X <- matrix(runif(200), nrow = 100)
  Y <- matrix(runif(200), nrow = 100)
  deco_features <- DecoELA(X = X, Y = Y, aggregate = TRUE, H = 5, scalar_func = "weightedsum", set_name = "disp")
  expect_is(deco_features, "list")
  expect_true(all(grepl("^deco\\.", names(deco_features))))
})

test_that("DecoELA correctly calculates features (nbc)", {
  X <- matrix(runif(200), nrow = 100)
  Y <- matrix(runif(200), nrow = 100)
  deco_features <- DecoELA(X = X, Y = Y, aggregate = TRUE, H = 5, scalar_func = "weightedsum", set_name = "nbc")
  expect_is(deco_features, "list")
  expect_true(all(grepl("^deco\\.", names(deco_features))))
})

test_that("DecoELA correctly calculates features (pca)", {
  X <- matrix(runif(200), nrow = 100)
  Y <- matrix(runif(200), nrow = 100)
  deco_features <- DecoELA(X = X, Y = Y, aggregate = TRUE, H = 5, scalar_func = "weightedsum", set_name = "pca")
  expect_is(deco_features, "list")
  expect_true(all(grepl("^deco\\.", names(deco_features))))
})

test_that("DecoELA correctly calculates features (ic)", {
  X <- matrix(runif(200), nrow = 100)
  Y <- matrix(runif(200), nrow = 100)
  deco_features <- DecoELA(X = X, Y = Y, aggregate = TRUE, H = 5, scalar_func = "weightedsum", set_name = "ic")
  expect_is(deco_features, "list")
  expect_true(all(grepl("^deco\\.", names(deco_features))))
})

test_that("DecoELA correctly calculates features (fdc)", {
  X <- matrix(runif(200), nrow = 100)
  Y <- matrix(runif(200), nrow = 100)
  deco_features <- DecoELA(X = X, Y = Y, aggregate = TRUE, H = 5, scalar_func = "weightedsum", set_name = "fdc")
  expect_is(deco_features, "list")
  expect_true(all(grepl("^deco\\.", names(deco_features))))
})

# Test for DecoELA with errors
test_that("DecoELA handles errors correctly", {
  X <- matrix(runif(20), nrow = 10)
  Y <- matrix(runif(20), nrow = 10)
  expect_error(DecoELA(X = "X", Y = Y, aggregate = TRUE, H = 5, scalar_func = "weightedsum", set_name = "ela_meta"))
  expect_error(DecoELA(X = X, Y = "Y", aggregate = TRUE, H = 5, scalar_func = "weightedsum", set_name = "ela_meta"))
  expect_error(DecoELA(X = X, Y = Y, aggregate = TRUE, H = 5, scalar_func = "invalid_scalar_func", set_name = "ela_meta"))
})

# Test case for empty matrix
test_that("DecoELA handles empty matrices", {
  X_empty <- matrix(numeric(0), nrow = 0, ncol = 0)
  Y_empty <- matrix(numeric(0), nrow = 0, ncol = 0)
  expect_error(DecoELA(X = X_empty, Y = Y_empty, aggregate = TRUE, H = 5, scalar_func = "weightedsum", set_name = "ela_meta"))
})

# Test case for mismatched rows in X and Y
test_that("DecoELA throws error for mismatched rows", {
  X <- matrix(runif(20), nrow = 10)
  Y <- matrix(runif(30), nrow = 15)
  expect_error(DecoELA(X = X, Y = Y, aggregate = TRUE, H = 5, scalar_func = "weightedsum", set_name = "ela_meta"))
})

# Test case for scalar_func being NULL
test_that("DecoELA throws error when scalar_func is NULL", {
  X <- matrix(runif(20), nrow = 10)
  Y <- matrix(runif(20), nrow = 10)

  expect_error(DecoELA(X = NULL, Y = Y, aggregate = TRUE, H = 5, scalar_func = NULL, set_name = "ela_meta"))
  expect_error(DecoELA(X = X, Y = NULL, aggregate = TRUE, H = 5, scalar_func = "weightedsum", set_name = "ela_meta"))
  expect_error(DecoELA(X = X, Y = Y, aggregate = TRUE, H = 5, scalar_func = "weightedsum", set_name = NULL))
})

