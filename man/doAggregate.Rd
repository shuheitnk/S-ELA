\name{doAggregate}
\alias{doAggregate}
\title{doAggregate: Aggregate Sub-Features}
\description{
  This function aggregates sub-feature sets by calculating descriptive statistics such as mean, min, max, and standard deviation in the decomposition-based approach.
}
\usage{
  doAggregate(sub_features)
}
\arguments{
  \item{sub_features}{A list of sub-features.}
}
\value{
  A list containing the aggregated features with statistics like mean, min, max, and standard deviation for each metric.
}
\examples{
  # Example usage of doAggregate
  sub_features <- list(list(metric1 = c(1, 2, 3), metric2 = c(4, 5, 6)),
                       list(metric1 = c(7, 8, 9), metric2 = c(10, 11, 12)))
  aggregated <- doAggregate(sub_features)
}
