library(testthat)
library(ScalarELA)

# Test for DecoELA function
test_that("DecoELA correctly calculates features", {
  # Create input data
  X <- matrix(runif(200), nrow = 100)
  Y <- matrix(runif(200), nrow = 100)

  # Check if DecoELA function runs without errors
  deco_features <- DecoELA(X = X, Y = Y, aggregate = TRUE, H = 5, make_vec = "sld", scalar_func = "weightedsum", set_name = "ela_meta")

  # Verify the output is a list
  expect_is(deco_features, "list", info = "The output of DecoELA should be a list.")

  # Check if the required fields are present in the output (example fields)
  expect_true("features" %in% names(deco_features), info = "The output of DecoELA should contain a 'features' field.")
  expect_true("meta_features" %in% names(deco_features), info = "The output of DecoELA should contain a 'meta_features' field.")
})

# Test for DecoELA error handling with invalid input
test_that("DecoELA handles errors correctly", {
  X <- matrix(runif(200), nrow = 100)
  Y <- matrix(runif(200), nrow = 100)

  # Test if an error is thrown when an invalid value is passed for H
  expect_error(DecoELA(X = X, Y = Y, aggregate = TRUE, H = "invalid", make_vec = "sld", scalar_func = "weightedsum", set_name = "ela_meta"),
               "H must be an integer", info = "DecoELA should throw an error if H is not a valid integer.")

  # Test if an error is thrown when an invalid scalar_func is provided
  expect_error(DecoELA(X = X, Y = Y, aggregate = TRUE, H = 5, make_vec = "sld", scalar_func = "invalid_func", set_name = "ela_meta"),
               "Invalid scalarization function", info = "DecoELA should throw an error for an invalid scalar_func.")
})
