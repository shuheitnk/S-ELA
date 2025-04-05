# Function to compute the Fitness-Distance Correlation (FDC)
computefdc <- function(X, Y) {
  flacco::measureTime(expression({

    # Extract the optimal solutions from the given set Y
    min_y <- min(Y)  # Find the minimum value in Y
    optima <- X[Y == min_y, , drop = FALSE]  # Extract all X values where Y is the minimum

    # Calculate the distance from each sample in X to the nearest optimal solution
    distances <- apply(X, 1, function(xi) {
      min(apply(optima, 1, function(o) sqrt(sum((xi - o)^2))))  # Find minimum distance
    })

    # Compute the Fitness-Distance Correlation (FDC)
    # Pearson correlation between objective values (Y) and distances to the nearest optimal solution
    fdc = cor(Y, distances, method = "pearson")

    list(fdc = fdc)
  }), "fitness_distance")
}
