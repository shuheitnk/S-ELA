# Function for decomposition-based approach
DecoELA = function(X, Y, H = 50, aggregate = TRUE, set_name = NULL){

  if (!requireNamespace("flacco", quietly = TRUE)) {
    stop("The 'flacco' package is required but not installed.")
  }
  
  # Ensure the number of solutions and fitness values match
  n.solutions = nrow(X)
  n.fitness = nrow(Y)
  
  if (n.solutions != n.fitness) {
    stop("Error: The number of solutions and the number of fitness values must be equal.")
  } else {
    
    # Weight vector for decomposition
    weight_vector <- decomposition_sld(list(name = "sld", H = H, .nobj = n.fitness))
    n.weight <- length(weight_vector)
    
    # Initialize result list
    sub_features <- vector("list", n.weight)
    
    # Loop over each weight to calculate the new objective values
    for (i in 1:n.weight) {
      w <- weight_vector[i]
      
      # Calculate new objective values based on the weight
      Fn <- sapply(1:nrow(Y), function(i) {
        w[i] * Y[i, j]  
      })
      
      feat_object <- flacco::createFeatureObject(X = X, y = Fn)
      
      # Calculate feature set for each weight
      sub_features[[i]] <- flacco::calculateFeatureSet(feat_object, set = set_name)
    }
    
    # Aggregate the results if required
    if (aggregate) {
      deco <- doAggregate(sub_features)
    } else {
      deco <- sub_features
    }
    
    return(deco)
  }
}

# Function for NDS (non-dominated sorting)-based approach
DomiELA = function(X, Y, set_name = NULL){

  if (!requireNamespace("flacco", quietly = TRUE)) {
    stop("The 'flacco' package is required but not installed.")
  }

  if (!requireNamespace("ecr", quietly = TRUE)) {
    stop("The 'ecr' package is required but not installed.")
  }
  
  # Ensure the number of solutions and fitness values match
  n.solutions = nrow(X)
  n.fitness = nrow(Y)
  
  if (n.solutions != n.fitness) {
    stop("Error: The number of solutions and the number of fitness values must be equal.")
  } else {
    
    # Perform non-dominated sorting and calculate feature set
    ranks <- ecr:::doNondominatedSortingR(t(Y))$ranks
    ranks <- as.numeric(ranks)
    
    feat_object <- flacco::createFeatureObject(X = X, y = ranks)
    
    # Calculate the dominance-based feature set
    domi <- flacco::calculateFeatureSet(feat_object, set = set_name)
    
    return(domi)
  }
}

  
  
  
