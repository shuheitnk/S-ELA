results_arrange <- function(results) {
  metrics <- names(results[[1]])
  
  # Initialize a list to store the results
  stat_results <- list()
  
  # Calculate statistics for each metric
  for (i in seq_along(metrics)) {
    values <- sapply(results[1:2], function(x) x[[metrics[i]]])
    
    # Calculate statistics for each metric and dynamically name them to store in the list
    stat_results[[paste0("deco.", metrics[i], ".mean")]] <- mean(values)
    stat_results[[paste0("deco.", metrics[i], ".min")]] <- min(values)
    stat_results[[paste0("deco.", metrics[i], ".max")]] <- max(values)
    stat_results[[paste0("deco.", metrics[i], ".sd")]] <- sd(values)
  }
  
  # Return the results as a list
  return(stat_results)
}



calculate_stats <- function(data) {
  
  # NAを除去
  data <- na.omit(data)
  n <- length(data)
  
  # 入力データのチェック
  if (n == 0) {
    stop("Error: adj_r2_data must not be empty or all NA.")
  }
  
  
  # 一度の走査で必要な統計量を計算
  min_val <- Inf
  max_val <- -Inf
  sum_val <- 0
  sum_sq_val <- 0
  
  # ソートされたデータを保持
  sorted_data <- numeric(n)
  
  for (i in seq_len(n)) {
    value <- data[i]
    sorted_data[i] <- value
    
    if (value < min_val) min_val <- value
    if (value > max_val) max_val <- value
    sum_val <- sum_val + value
    sum_sq_val <- sum_sq_val + value^2
  }
  
  # ソート
  sorted_data <- sort(sorted_data)
  
  mean_val <- sum_val / n
  sd_val <- sqrt((sum_sq_val / n) - (mean_val^2))
  median_val <- median(sorted_data)
  
  # 結果をリストにまとめる
  result <- list(
    min = min_val,
    mean = mean_val,
    max = max_val,
    sd = sd_val,
    med_mean_diff = median_val - mean_val
  )
  
  return(result)
}


minmax <- function(data) {
  # NAを除去
  data <- na.omit(data)
  n <- length(data)
  # データのチェック
  if (n == 0) {
    stop("Error: adj_r2_data must not be empty or all NA.")
  }
  
  # 一度の走査で平均と標準偏差を計算
  
  sum_val <- 0
  sum_sq_val <- 0
  
  for (value in data) {
    sum_val <- sum_val + value
    sum_sq_val <- sum_sq_val + value^2
  }
  
  mean_val <- sum_val / n
  std_dev <- sqrt((sum_sq_val / n) - (mean_val^2))
  obj_scale_cv = std_dev / mean_val
  
  return(obj_scale_cv)
}


calculateDistributionFeatures = function(y) {
  
  measureTime(expression({
    list(ela_distr.skewness = e1071::skewness(y),
         ela_distr.kurtosis = e1071::kurtosis(y))
  }), "ela_distr")
}
