library(testthat)
library(ScalarELA)

# Test for DomiELA function
test_that("DomiELA correctly calculates features", {
  # Create input data
  X <- matrix(runif(200), nrow = 100)
  Y <- matrix(runif(200), nrow = 100)

  # Check if DomiELA function runs without errors
  domi_features <- DomiELA(X = X, Y = Y, aggregate = TRUE, make_vec = "sld", scalar_func = "weightedsum", set_name = "ela_meta")

  # Verify the output is a list
  expect_is(domi_features, "list", info = "The output of DomiELA should be a list.")

  # Check if the required fields are present in the output (example fields)
  expect_true("features" %in% names(domi_features), info = "The output of DomiELA should contain a 'features' field.")
  expect_true("meta_features" %in% names(domi_features), info = "The output of DomiELA should contain a 'meta_features' field.")
})

# Test for DomiELA error handling with invalid input
test_that("DomiELA handles errors correctly", {
  X <- matrix(runif(200), nrow = 100)
  Y <- matrix(runif(200), nrow = 100)

  # Test if an error is thrown when an invalid value is passed for aggregate
  expect_error(DomiELA(X = X, Y = Y, aggregate = "invalid", make_vec = "sld", scalar_func = "weightedsum", set_name = "ela_meta"),
               "aggregate should be a logical value", info = "DomiELA should throw an error if aggregate is not a logical value.")

  # Test if an error is thrown when an invalid scalar_func is provided
  expect_error(DomiELA(X = X, Y = Y, aggregate = TRUE, make_vec = "sld", scalar_func = "invalid_func", set_name = "ela_meta"),
               "Invalid scalarization function", info = "DomiELA should throw an error for an invalid scalar_func.")
})
