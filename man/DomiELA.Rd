\name{DomiELA}
\alias{DomiELA}
\title{DomiELA: Feature Calculation for NDS-based Approach}
\description{
  This function calculates ELA feature sets based on the non-dominated sorting (NDS) in multi-objective optimization.
}
\usage{
  DomiELA(X, Y, set_name)
}
\arguments{
  \item{X}{A matrix or data frame containing the sampled solutions (rows are solutions, columns are variables).}
  \item{Y}{A matrix or data frame containing the fitness values of the sampled solutions (rows are solutions, columns are objectives).}
  \item{set_name}{A string specifying ELA feature set to be used.}
}
\value{
  A list containing the S-ELA features calculated from the NDS-based approach.
}
\examples{
  # Example usage of DomiELA
  X <- matrix(runif(200), nrow = 100)
  Y <- matrix(runif(200), nrow = 100)
  domi_features <- DomiELA(X = X,
                           Y = Y,
                           set_name = "ela_meta")

}
