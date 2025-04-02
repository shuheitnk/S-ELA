\name{weightedsum}
\alias{weightedsum}
\title{weightedsum: Weighted Sum Scalarization Function}
\description{
  This function calculates the scalarized objective values using the weighted sum method.
}
\usage{
  weightedsum(Y, w)
}
\arguments{
  \item{Y}{A matrix or data frame containing the fitness values of the sampled solutions (rows are solutions, columns are objectives).}
  \item{w}{A weight vector containing the weights for convex combination. The number of weights must match the number of objectives.}
}
\value{
  A vector of scalarized objective values for each solution.
}
\examples{
  # Example usage of weightedsum
  Y <- matrix(runif(20), nrow = 10)
  w <- c(0.3, 0.7)
  scalarized_values <- weightedsum(Y, w)
}
