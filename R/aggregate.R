doAggregate <- function(sub_features) {
  metrics <- names(sub_features[[1]])
  aggregated_features <- list()
  
  for (metric in metrics) {
    values <- sapply(sub_features[1:2], function(x) x[[metric]])
    
    aggregated_features[[paste0("deco.", metric, ".mean")]] <- mean(values)
    aggregated_features[[paste0("deco.", metric, ".min")]] <- min(values)
    aggregated_features[[paste0("deco.", metric, ".max")]] <- max(values)
    aggregated_features[[paste0("deco.", metric, ".sd")]] <- sd(values)
  }
  
  return(aggregated_features)
}
