library(testthat)
library(ScalarELA)
library(flacco)
if (!requireNamespace("RANN", quietly = TRUE)) {
  install.packages("RANN")
}
library(RANN)

test_check("ScalarELA")
