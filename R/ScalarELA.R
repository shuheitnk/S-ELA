# Function for decomposition-based approach
DecoELA = function(X, Y, H = 50, aggregate = TRUE, make_vec = "sld", scalar_func = "weightedsum", set_name = NULL) {
  
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
    weight_vector <- MOEADr::decomposition_sld(list(name = make_vec, H = H, .nobj = ncol(Y)))
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

# Function for NDS (non-dominated sorting)-based approach
DomiELA = function(X, Y, set_name = NULL) {
  
  # Check if required packages are installed
  if (!requireNamespace("flacco", quietly = TRUE)) {
    stop("The 'flacco' package is required but not installed. Please install it using install.packages('flacco').")
  }
  
  if (!requireNamespace("ecr", quietly = TRUE)) {
    stop("The 'ecr' package is required but not installed. Please install it using install.packages('ecr').")
  }
  
  # Ensure the number of solutions and fitness values match
  n.solutions = nrow(X)
  n.fitness = nrow(Y)
  
  if (n.solutions != n.fitness) {
    stop("Error: The number of solutions (", n.solutions, ") and fitness values (", n.fitness, ") must be equal.")
  } else {
    
    # Perform non-dominated sorting on the fitness values
    ranks <- ecr:::doNondominatedSortingR(t(Y))$ranks
    R <- as.numeric(ranks)
    
    # Create feature object for the non-dominated solutions
    feat_object <- flacco::createFeatureObject(X = X, y = R)
    
    # Calculate the dominance-based feature set
    domi <- flacco::calculateFeatureSet(feat_object, set = set_name)
    
    # Process feature set into a list with modified names
    metrics <- names(domi)
    domi_features <- list()
    
    for (metric in metrics) {
      values <- domi[[metric]]
      domi_features[[paste0("domi.", gsub("^ela_", "", metric))]] <- values
    }
    
    return(domi_features)
  }
}
