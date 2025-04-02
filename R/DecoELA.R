# Function for decomposition-based approach
DecoELA = function(X, Y, H, aggregate = TRUE, scalar_func = "weightedsum", set_name) {

  # Check if required packages are installed
  if (!requireNamespace("flacco", quietly = TRUE)) {
    stop("The 'flacco' package is required but not installed. Please install it using install.packages('flacco').")
  }

  if (!requireNamespace("MOEADr", quietly = TRUE)) {
    stop("The 'MOEADr' package is required but not installed. Please install it using install.packages('MOEADr').")
  }

  # Select scalarization function based on user input
  if (scalar_func == "weightedsum") {
    Scfunc = weightedsum
  } else if (scalar_func == "tchebycheff") {
    Scfunc = tchebycheff
  } else {
    stop("Error: The specified scalarization function ('", scalar_func, "') is not implemented.")
  }

  # Ensure the number of solutions and fitness values match
  n.solutions = nrow(X)
  n.fitness = nrow(Y)

  if (n.solutions != n.fitness) {
    stop("Error: The number of solutions (", n.solutions, ") and fitness values (", n.fitness, ") must be equal.")
  } else {

    # Generate the weight vector for decomposition
    weight_vector <- MOEADr::decomposition_sld(list(name = "sld", H = H, .nobj = ncol(Y)))
    n.weight <- nrow(weight_vector)

    # Initialize the result list for features
    sub_features <- vector("list", n.weight)

    # Loop over each weight vector to calculate the new objective values and feature sets
    for (i in 1:n.weight) {
      w <- weight_vector[i,]

      # Compute the scalarized objective values based on the weight vector
      G <- Scfunc(Y, w)

      # Create feature object for the current set of solutions
      feat_object <- flacco::createFeatureObject(X = X, y = G)

      # Calculate the feature set for the current weight vector
      sub_features[[i]] <- flacco::calculateFeatureSet(feat_object, set = set_name)
    }

    # Aggregate the feature sets if required
    if (aggregate) {
      deco_features <- doAggregate(sub_features)
    } else {
      deco_features <- sub_features
    }

    return(deco_features)
  }
}
