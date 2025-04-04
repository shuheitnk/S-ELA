test_that("DomiELA correctly calculates features (ela_meta)", {
  X <- matrix(runif(200), nrow = 100)
  Y <- matrix(runif(200), nrow = 100)
  domi_features = DomiELA(X, Y, set_name = "ela_meta")
  expect_is(domi_features, "list")
  expect_true(all(grepl("^domi\\.", names(domi_features))))
})


test_that("DomiELA correctly calculates features (ela_distr)", {
  X <- matrix(runif(200), nrow = 100)
  Y <- matrix(runif(200), nrow = 100)
  domi_features = DomiELA(X, Y, set_name = "ela_distr")
  expect_is(domi_features, "list")
  expect_true(all(grepl("^domi\\.", names(domi_features))))
})

test_that("DomiELA correctly calculates features (ela_meta, basic)", {

  X <- matrix(runif(200), nrow = 100)
  Y <- matrix(runif(200), nrow = 100)
  domi_features = DomiELA(X, Y, set_name = "basic")
  expect_is(domi_features, "list")
  expect_true(all(grepl("^domi\\.", names(domi_features))))
})

test_that("DomiELA correctly calculates features (disp)", {

  X <- matrix(runif(200), nrow = 100)
  Y <- matrix(runif(200), nrow = 100)
  domi_features = DomiELA(X, Y, set_name = "disp")
  expect_is(domi_features, "list")
  expect_true(all(grepl("^domi\\.", names(domi_features))))
})


test_that("DomiELA correctly calculates features (nbc)", {

  X <- matrix(runif(200), nrow = 100)
  Y <- matrix(runif(200), nrow = 100)
  domi_features = DomiELA(X, Y, set_name = "nbc")
  expect_is(domi_features, "list")
  expect_true(all(grepl("^domi\\.", names(domi_features))))
})

test_that("DomiELA correctly calculates features (pca)", {

  X <- matrix(runif(200), nrow = 100)
  Y <- matrix(runif(200), nrow = 100)
  domi_features = DomiELA(X, Y, set_name = "pca")
  expect_is(domi_features, "list")
  expect_true(all(grepl("^domi\\.", names(domi_features))))
})


test_that("DomiELA correctly calculates features (ic)", {

  X <- matrix(runif(200), nrow = 100)
  Y <- matrix(runif(200), nrow = 100)
  domi_features = DomiELA(X, Y, set_name = "ic")

  expect_is(domi_features, "list")
  expect_true(all(grepl("^domi\\.", names(domi_features))))
})




test_that("DomiELA handles errors correctly", {
  X <- matrix(runif(20), nrow = 10)
  Y <- matrix(runif(20), nrow = 10)
  expect_error(DomiELA(X = "X", Y = Y, set_name = "ela_meta"))
  expect_error(DomiELA(X = X, Y = "Y", set_name = "ela_meta"))
  expect_error(DomiELA(X = X, Y = Y, set_name = "invalid_set"))
})

# Test case for empty matrix
test_that("DomiELA handles empty matrices", {
  X_empty <- matrix(numeric(0), nrow = 0, ncol = 0)
  Y_empty <- matrix(numeric(0), nrow = 0, ncol = 0)
  expect_error(DomiELA(X = X_empty, Y = Y_empty, set_name = "ela_meta"))

})

# Test case for different number of rows in X and Y
test_that("DomiELA throws error for mismatched rows", {
  X <- matrix(runif(20), nrow = 10)
  Y <- matrix(runif(30), nrow = 15)
  expect_error(DomiELA(X = X, Y = Y, set_name = "ela_meta"))
})

# Test case for scalar_func being NULL
test_that("DomiELA throws error when scalar_func is NULL", {
  X <- matrix(runif(20), nrow = 10)
  Y <- matrix(runif(20), nrow = 10)

  expect_error(DomiELA(X = NULL, Y = Y, set_name = "ela_meta"))
  expect_error(DomiELA(X = X, Y = NULL, set_name = "ela_meta"))
  expect_error(DomiELA(X = X, Y = Y, set_name = NULL))
})
