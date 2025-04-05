

# Function for NDS (non-dominated sorting)-based approach
DomiELA = function(X, Y, normalize_X = TRUE, normalize_R = TRUE, set_name) {

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
    ranks <- ecr::doNondominatedSorting(t(Y))$ranks
    R <- as.numeric(ranks)

    if (normalize_X == TRUE){
      X <- NormalizeColumns(X)
    }

    if (normalize_R == TRUE){
      R <- NormalizeColumns(matrix(R, ncol = 1))
    }

    # Create feature object for the non-dominated solutions
    feat_object <- flacco::createFeatureObject(X = X, y = R)

    # Calculate the dominance-based feature set
    if (set_name == "ela_distr") {
      domi <- domiDistr(Y)
    }else if (set_name == "fdc") {
      domi <- computefdc(X = X, Y = R)
    }else{
      domi <- flacco::calculateFeatureSet(feat_object, set = set_name)
    }



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
