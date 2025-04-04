\name{domiDistr}
\alias{domiDistr}
\title{domiDistr: ela_distr in NDS-based approach}
\description{
  This function calculates the skewness and kurtosis of a given vector of values
  in the context of the non-dominated sorting (NDS) approach used in multi-objective optimization.
  It was developed to handle scenarios where errors were encountered when calculating the number of peaks during ranking.
}
\usage{
  domiDistr(y, skewness.type = 3L, kurtosis.type = 3L)
}
\arguments{
  \item{y}{A numeric vector}
  \item{skewness.type}{The type of skewness to compute (default is 3L).}
  \item{kurtosis.type}{The type of kurtosis to compute (default is 3L).}

}
\value{
  A list of ela_distr features.
}
\examples{
  # Example usage of domiDistr

  Y <- matrix(runif(100), nrow = 100)
  distr <- domiDistr(y = Y)
}

