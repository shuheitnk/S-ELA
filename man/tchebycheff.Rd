\name{tchebycheff}
\alias{tchebycheff}
\title{tchebycheff: Tchebycheff Scalarization Function}
\description{
  This function calculates the scalarized objective values using the Tchebycheff method.
}
\usage{
  tchebycheff(Y, w)
}
\arguments{
  \item{Y}{A matrix or data frame containing the fitness values (rows are solutions, columns are objectives).}
  \item{w}{A weight vector containing the weights for the Tchebycheff method. The number of weights must match the number of objectives.}
}
\value{
  A vector of scalarized objective values for each solution.
}
\examples{
  # Example usage of tchebycheff
  Y <- matrix(runif(20), nrow = 10)
  w <- c(0.5, 0.5)
  scalarized_values <- tchebycheff(Y, w)
}
