# Function to generate sample points and evaluate the objective function
makeObject <- function(fn, n.sample) {
  
  # Check for required libraries and stop if not installed
  if (!requireNamespace("smoof", quietly = TRUE)) {
    stop("The 'smoof' package is required but not installed.")
  }
  if (!requireNamespace("lhs", quietly = TRUE)) {
    stop("The 'lhs' package is required but not installed.")
  }
  
  # Get the lower and upper bounds of the function's box constraints
  fn.lower <- smoof::getLowerBoxConstraints(fn)
  fn.upper <- smoof::getUpperBoxConstraints(fn)
  
  # Number of decision variables
  d <- length(fn.lower)
  
  # Ensure the function has at least one decision variable
  if (d < 1) {
    stop("The function must have at least one decision variable.")
  }
  
  # Generate initial Latin Hypercube Samples (LHS)
  samples <- lhs::improvedLHS(n.sample, d)
  
  # Scale samples to the defined bounds of the function
  samples <- sweep(samples, 2, fn.lower, FUN = "*") + fn.lower
  
  # Evaluate the objective function at the generated samples
  y <- t(apply(samples, 1, fn))
  
  # If the output is a vector, convert it to a matrix
  if (is.vector(y)) {
    y <- matrix(y, ncol = 1)
  }
  
  # Number of objectives (columns of the output)
  m <- ncol(y)
  
  # Scale the objectives to the [0, 1] range
  y_min <- apply(y, 2, min)
  y_max <- apply(y, 2, max)
  obj_scale <- y_max - y_min
  y <- sweep(y, 2, y_min, FUN = "-") / obj_scale
  
  # Scale the samples to the [0, 1] range
  samples <- sweep(samples, 2, fn.lower, FUN = "-") / (fn.upper - fn.lower)
  
  # Return the scaled samples, objectives, and scaling factors
  return(list(samples = samples, y = y, obj_scale = obj_scale))
}

