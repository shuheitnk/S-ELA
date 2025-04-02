# Function to aggregate sub-features (e.g., mean, min, max, standard deviation)
doAggregate <- function(sub_features) {

  # Ensure that sub_features is not empty and has at least one element
  if (length(sub_features) == 0 || !is.list(sub_features[[1]])) {
    stop("Error: 'sub_features' must be a non-empty list with valid sub-feature structures.")
  }

  # Extract the names of the metrics from the first element of sub_features
  metrics <- names(sub_features[[1]])

  # Initialize an empty list to store aggregated features
  aggregated_features <- list()

  # Loop over each metric to calculate the aggregation statistics
  for (metric in metrics) {
    # Extract the values for the current metric from all sub-feature sets
    values <- sapply(sub_features, function(x) x[[metric]])

    # Calculate mean, min, max, and standard deviation for the metric
    aggregated_features[[paste0("deco.", gsub("^ela_", "", metric), ".mean")]] <- mean(values)
    aggregated_features[[paste0("deco.", gsub("^ela_", "", metric), ".min")]] <- min(values)
    aggregated_features[[paste0("deco.", gsub("^ela_", "", metric), ".max")]] <- max(values)
    aggregated_features[[paste0("deco.", gsub("^ela_", "", metric), ".sd")]] <- sd(values)
  }

  # Return the aggregated features
  return(aggregated_features)
}
