domiDistr = function(y, skewness.type = 3L, kurtosis.type = 3L){

  # A version of the function has been developed that calculates only the skewness and kurtosis,
  # as errors were encountered when attempting to compute the number of peaks during ranking.
  # This modification resolves the issue by focusing exclusively on the distributional properties
  # (skewness and kurtosis), without including the calculation of peaks.

  flacco::measureTime(expression({
    list(ela_distr.skewness = e1071::skewness(y, type = skewness.type),
         ela_distr.kurtosis = e1071::kurtosis(y, type = kurtosis.type))
  }), "ela_distr")

}
