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
  
  if (length(fn.lower) < 2 || length(fn.upper) < 2) {
    stop("The function must have at least two dimensions for x1 and x2.")
  }
  
  # Set ranges for x1 and x2
  x1_range <- c(fn.lower[1], fn.upper[1])  # Range for x1
  x2_range <- c(fn.lower[2], fn.upper[2])  # Range for x2
  
  # Generate initial samples
  samples <- lhs::improvedLHS(n.sample, 2)
  
  # Scale samples to the defined ranges
  samples <- (samples * c(diff(x1_range), diff(x2_range))) + c(x1_range[1], x2_range[1])
  
  # Apply the function to samples
  y <- t(apply(samples, 1, fn))  # Transpose to ensure correct structure
  
  # Split y into y1 and y2
  y1 <- y[, 1]
  y2 <- y[, 2]
  
  # Calculate min and max for scaling
  y1_min <- min(y1)
  y1_max <- max(y1)
  y2_min <- min(y2)
  y2_max <- max(y2)
  
  # Scale y1 and y2
  y1 <- (y1 - y1_min) / (y1_max - y1_min)
  y2 <- (y2 - y2_min) / (y2_max - y2_min)
  
  # Calculate objective scaling
  obj_scale <- c((y1_max - y1_min), (y2_max - y2_min))
  
  # Update y with scaled values
  y[, 1] <- y1
  y[, 2] <- y2
  
  # Scale samples to [0, 1]
  samples[, 1] <- (samples[, 1] - x1_range[1]) / diff(x1_range)
  samples[, 2] <- (samples[, 2] - x2_range[1]) / diff(x2_range)
  
  # Return samples and y
  return(list(samples = samples, y = y, obj_scale = obj_scale))
}


