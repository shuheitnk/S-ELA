makeobject <- function(fn, n.sample) {
  # Load required libraries
  if (!requireNamespace("smoof", quietly = TRUE)) {
    stop("The 'smoof' package is required but not installed.")
  }
  if (!requireNamespace("lhs", quietly = TRUE)) {
    stop("The 'lhs' package is required but not installed.")
  }
  
  # Get the box constraints of the function
  fn.lower <- smoof::getLowerBoxConstraints(fn)
  fn.upper <- smoof::getUpperBoxConstraints(fn)
  
  d <- length(fn.lower)  # Number of decision variables
  if (d < 1) {
    stop("The function must have at least one decision variable.")
  }
  
  # Generate initial samples
  samples <- lhs::improvedLHS(n.sample, d)
  
  # Scale samples to the defined ranges
  samples <- sweep(samples, 2, fn.lower, FUN = "*") + fn.lower
  
  # Evaluate the function
  y <- t(apply(samples, 1, fn))
  
  # Ensure y is a matrix
  if (is.vector(y)) {
    y <- matrix(y, ncol = 1)
  }
  
  m <- ncol(y)  # Number of objectives
  
  # Scale objectives
  y_min <- apply(y, 2, min)
  y_max <- apply(y, 2, max)
  obj_scale <- y_max - y_min
  y <- sweep(y, 2, y_min, FUN = "-") / obj_scale
  
  # Scale samples to [0, 1]
  samples <- sweep(samples, 2, fn.lower, FUN = "-") / (fn.upper - fn.lower)
  
  # Return results
  return(list(samples = samples, y = y, obj_scale = obj_scale))
}

