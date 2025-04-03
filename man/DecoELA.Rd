\name{DecoELA}
\alias{DecoELA}
\title{DecoELA: Feature Calculation for Decomposition-based Approach}
\description{
  This function calculates ELA feature sets based on the decomposition approach in multi-objective optimization.
}
\usage{
  DecoELA(X, Y, H, aggregate = TRUE, scalar_func = "weightedsum", set_name)
}
\arguments{
  \item{X}{A matrix or data frame containing the sampled solutions (rows are solutions, columns are variables).}
  \item{Y}{A matrix or data frame containing the fitness values of the sampled solution (rows are solutions, columns are objectives).}
  \item{H}{An integer specifying the number of divisions for decomposition. }
  \item{aggregate}{A logical value indicating whether to aggregate sub-feature sets. Default is TRUE.}
  \item{scalar_func}{A string specifying the scalarization function to use. Options are "weightedsum" or "tchebycheff". Default is "weightedsum".}
  \item{set_name}{A string specifying ELA feature set to be used. }
}
\value{
  A list containing the S-ELA features calculated from the decomposition-based approach.
}
\examples{
  # Example usage of DecoELA
  X <- matrix(runif(200), nrow = 100)
  Y <- matrix(runif(200), nrow = 100)
  deco_features <- DecoELA(X = X,
                           Y = Y,
                           aggregate = TRUE,
                           H = 5,
                           scalar_func = "weightedsum",
                           set_name = "ela_meta")
}
