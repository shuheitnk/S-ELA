\name{computefdc}
\alias{computefdc}
\title{
  A function to compute the Fitness-Distance Correlation (FDC)
}
\description{
  This function calculates the Fitness-Distance Correlation (FDC), which quantifies the relationship between the objective values of solutions and their distance to the nearest optimal solution.
}
\usage{
  computefdc(X, Y)
}
\arguments{
  \item{X}{
    A numeric matrix of sampled solutions.
  }
  \item{Y}{
    A numeric vector of the objective function values for each solution in \code{X}.
  }
}
\details{
  The Fitness-Distance Correlation (FDC) is a measure that helps assess how well the fitness (objective function value) correlates with the distance to the global optimum. The function works as follows:

  1. Identifies the set of optimal solutions, defined as the solutions that have the minimum objective value in \code{Y}.
  2. For each solution in \code{X}, computes the Euclidean distance to the closest optimal solution.
  3. Computes the Pearson correlation between the fitness values \code{Y} and the calculated distances.
}
\value{
  A numeric value representing the Pearson correlation coefficient (r) between the objective values \code{Y} and the distances from each solution to the nearest optimal solution.
}
\examples{
  ## Example of using compute_FDC function
  set.seed(123)
  X <- matrix(runif(100), nrow = 20, ncol = 5)  # 20 solutions, 5 variables each
  Y <- rowSums(X^2)  # Objective function: sum of squares (Sphere function)
  fdc_value <- computefdc(X, Y)
}
