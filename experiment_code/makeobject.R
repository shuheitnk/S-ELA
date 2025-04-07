makeobject <- function(fn, n.sample) {
  # Load required libraries
  if (!requireNamespace("smoof", quietly = TRUE)) {
    stop("The 'smoof' package is required but not installed.")
  }
  if (!requireNamespace("lhs", quietly = TRUE)) {
    stop("The 'lhs' package is required but not installed.")
  }
  
  # Retrieve lower and upper bounds of the decision space
  fn.lower <- smoof::getLowerBoxConstraints(fn)
  fn.upper <- smoof::getUpperBoxConstraints(fn)
  
  # Ensure the function is at least two-dimensional
  if (length(fn.lower) < 2 || length(fn.upper) < 2) {
    stop("The objective function must have at least two decision variables.")
  }
  
  # Define the ranges for the first two decision variables
  x1_range <- c(fn.lower[1], fn.upper[1])
  x2_range <- c(fn.lower[2], fn.upper[2])
  
  # Generate Latin Hypercube Samples in 2D
  samples <- lhs::improvedLHS(n.sample, 2)
  
  # Scale samples to the function's domain
  samples <- (samples * c(diff(x1_range), diff(x2_range))) + c(x1_range[1], x2_range[1])
  
  # Evaluate the objective function for each sample
  y <- t(apply(samples, 1, fn))  # Transpose to ensure each row corresponds to a solution
  
  # Separate objective values
  y1 <- y[, 1]
  y2 <- y[, 2]
  
  # Compute min and max for normalization
  y1_min <- min(y1)
  y1_max <- max(y1)
  y2_min <- min(y2)
  y2_max <- max(y2)
  
  # Normalize objective values to [0, 1]
  y1 <- (y1 - y1_min) / (y1_max - y1_min)
  y2 <- (y2 - y2_min) / (y2_max - y2_min)
  
  # Store original objective value ranges
  obj_scale <- c((y1_max - y1_min), (y2_max - y2_min))
  
  # Update the objective matrix with normalized values
  y[, 1] <- y1
  y[, 2] <- y2
  
  # Normalize decision variables to [0, 1] as well
  samples[, 1] <- (samples[, 1] - x1_range[1]) / diff(x1_range)
  samples[, 2] <- (samples[, 2] - x2_range[1]) / diff(x2_range)
  
  # Return the sample set, normalized objective values, and original objective ranges
  return(list(sample = samples, y = y, obj_scale = obj_scale))
}

