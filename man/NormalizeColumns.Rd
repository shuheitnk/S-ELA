\name{NormalizeColumns}
\alias{NormalizeColumns}
\title{
  Normalize Columns of a Matrix to the [0, 1] Range
}
\description{
  This function normalizes each column of the input matrix to the range [0, 1] by subtracting the minimum value of each column and dividing by the range (max - min) of the column.
}
\usage{
NormalizeColumns(Y)
}
\arguments{
  \item{Y}{
    A numeric matrix whose columns need to be normalized. Each column will be normalized individually.
  }
}
\details{
  The function applies min-max normalization to each column of the input matrix. For each column, the minimum value is subtracted, and the result is divided by the range (maximum value minus minimum value) of that column. This scales the values to the [0, 1] interval.
}
\value{
  A numeric matrix of the same dimensions as \code{Y}, but with each column normalized to the range [0, 1].
}
\examples{
# Example of normalizing a 5x4 matrix
Y <- matrix(runif(20), nrow = 5, ncol = 4)
Y_normalized <- NormalizeColumns(Y)

print(Y_normalized)
}
