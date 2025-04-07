# Aggregates statistics (mean, min, max, sd) for each feature across subproblems
rename_deco <- function(sub_features) {
  # Extract the feature names (assumed to be consistent across subproblems)
  metrics <- names(sub_features[[1]])
  
  # Initialize the output list
  aggregated_features <- list()
  
  # Loop over each feature and calculate summary statistics
  for (metric in metrics) {
    # Extract metric values from each subproblem
    values <- sapply(sub_features, function(x) x[[metric]])
    
    # Strip 'ela_' prefix (if any) and store statistics
    base_name <- gsub("^ela_", "", metric)
    aggregated_features[[paste0(base_name, ".mean")]] <- mean(values)
    aggregated_features[[paste0(base_name, ".min")]] <- min(values)
    aggregated_features[[paste0(base_name, ".max")]] <- max(values)
    aggregated_features[[paste0(base_name, ".sd")]]   <- stats::sd(values)
  }
  
  return(aggregated_features)
}

# Removes 'ela_' prefix from metric names without aggregation
rename_domi <- function(domi) {
  metrics <- names(domi)
  domi_features <- list()
  
  for (metric in metrics) {
    base_name <- gsub("^ela_", "", metric)
    domi_features[[base_name]] <- domi[[metric]]
  }
  
  return(domi_features)
}

# Computes basic statistics (min, mean, max, sd, median-mean difference)
calculate_stats <- function(data) {
  data <- na.omit(data)
  n <- length(data)
  
  if (n == 0) {
    stop("Error: 'data' must not be empty or consist only of NA values.")
  }
  
  # Initialize accumulators
  sum_val <- 0
  sum_sq_val <- 0
  min_val <- Inf
  max_val <- -Inf
  sorted_data <- numeric(n)
  
  # Compute aggregates
  for (i in seq_len(n)) {
    value <- data[i]
    sorted_data[i] <- value
    if (value < min_val) min_val <- value
    if (value > max_val) max_val <- value
    sum_val <- sum_val + value
    sum_sq_val <- sum_sq_val + value^2
  }
  
  # Compute final statistics
  sorted_data <- sort(sorted_data)
  mean_val <- sum_val / n
  sd_val <- sqrt((sum_sq_val / n) - mean_val^2)
  median_val <- median(sorted_data)
  
  return(list(
    min = min_val,
    mean = mean_val,
    max = max_val,
    sd = sd_val,
    med_mean_diff = median_val - mean_val
  ))
}

# Computes the coefficient of variation (CV = SD / mean)
minmax <- function(data) {
  data <- na.omit(data)
  n <- length(data)
  
  if (n == 0) {
    stop("Error: 'data' must not be empty or consist only of NA values.")
  }
  
  sum_val <- sum(data)
  sum_sq_val <- sum(data^2)
  
  mean_val <- sum_val / n
  std_dev <- sqrt((sum_sq_val / n) - mean_val^2)
  
  return(std_dev / mean_val)
}

# Computes skewness and kurtosis of a distribution using flacco's timing utility
calculateDistributionFeatures <- function(y) {
  flacco::measureTime(expression({
    list(
      ela_distr.skewness = e1071::skewness(y),
      ela_distr.kurtosis = e1071::kurtosis(y)
    )
  }), "ela_distr")
}

