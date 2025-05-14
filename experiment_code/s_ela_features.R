# Load required custom functions
source("Supportingfunctions.R")  #Several functions supporting calculation of S-ELA
source("makeobject.R")            # Function for generating input samples and objective values

# Set up experimental parameters
dimension <- 2
set.seed(100)
unique_seeds <- sample(1:1000, 1, replace = FALSE)  # Random seed for reproducibility
instance_ids <- 1:31  # BBOB instance IDs to evaluate

write_count <- 0  # Counter to track file writing

# Record the total time taken by the experiment
total_time <- system.time({
  
  instance_counter <- 0
  for (iid in instance_ids) {
    
    instance_counter <- instance_counter + 1
    
    # Define 55 bi-objective BBOB problems
    problems <- lapply(1:55, function(fid) {
      list(fn = smoof::makeBiObjBBOBFunction(dimensions = dimension, fid = fid, iid = iid))
    })
    
    for (problem in problems) {
      repetition_id <- 0
      
      for (seed in unique_seeds) {
        repetition_id <- repetition_id + 1
        
        feature_time <- system.time({
          
          set.seed(seed)
          
          obj_fn <- problem$fn
          instance_name <- smoof::getName(obj_fn)
          
          # Generate input samples and corresponding objective values
          obj_data <- makeobject(fn = obj_fn, n.sample = 200)
          X <- obj_data$sample
          Y <- obj_data$y
          obj_scale <- obj_data$obj_scale
          
          # Perform non-dominated sorting to assign ranks
          nds_ranks <- as.numeric(ecr:::doNondominatedSorting(t(Y))$ranks)
          num_samples <- length(nds_ranks)
          rank1_ratio <- sum(nds_ranks == 1) / num_samples
          
          # Generate weight vectors using MOEA/D decomposition
          weights <- MOEADr::decomposition_sld(list(name = "sld", H = 50, .nobj = 2))[, 1]
          num_weights <- length(weights)
          
          # Initialize feature object for NDS-based approach
          feat_obj_nds <- flacco::createFeatureObject(X = X, y = nds_ranks)
          
          # Prepare result lists
          results_meta <- results_distr <- results_disp <- results_nbc <- 
            results_ic <- results_fdc <- results_pca <- vector("list", num_weights)
          
          # Calculate decomposition-based features
          for (i in seq_along(weights)) {
            w <- weights[i]
            scalarized_Y <- sapply(1:nrow(Y), function(j) {
              w * Y[j, 1] + (1 - w) * Y[j, 2]
            })
            
            feat_obj <- flacco::createFeatureObject(X = X, y = scalarized_Y)
            
            results_meta[[i]]  <- flacco::calculateFeatureSet(feat_obj, set = "ela_meta")
            results_distr[[i]] <- flacco::calculateFeatureSet(feat_obj, set = "ela_distr")
            results_disp[[i]]  <- flacco::calculateFeatureSet(feat_obj, set = "disp")
            results_nbc[[i]]   <- flacco::calculateFeatureSet(feat_obj, set = "nbc")
            results_ic[[i]]    <- flacco::calculateFeatureSet(feat_obj, set = "ic")
            results_fdc[[i]]   <- flacco:::calculateFitnessDistanceFeatures(feat_obj) #fdc function in flacco
            results_pca[[i]]   <- flacco::calculateFeatureSet(feat_obj, set = "pca")
          }
          
          # Aggregate features into a single list
          feature_vector <- c(
            instance = instance_name,
            obj_scale = minmax(obj_scale),
            rank1_ratio = rank1_ratio,
            
            # Decomposition-based features
            deco = unlist(rename_deco(results_meta), use.names = TRUE),
            deco = unlist(rename_deco(results_distr)),
            deco = unlist(rename_deco(results_disp)),
            deco = unlist(rename_deco(results_nbc)),
            deco = unlist(rename_deco(results_ic)),
            deco = unlist(rename_deco(results_fdc)),
            deco = unlist(rename_deco(results_pca)),
            
            # NDS-based features (partially modified from flacco's original functions)
            domi = rename_domi(flacco::calculateFeatureSet(feat_obj_nds, set = "ela_meta")),
            domi = rename_domi(calculateDistributionFeatures(nds_ranks)),  # Modified function to handle rank values
            domi = rename_domi(flacco::calculateFeatureSet(feat_obj_nds, set = "disp")),
            domi = rename_domi(flacco::calculateFeatureSet(feat_obj_nds, set = "nbc")),
            domi = rename_domi(flacco::calculateFeatureSet(feat_obj_nds, set = "ic")),
            domi = rename_domi(flacco:::calculateFitnessDistanceFeatures(feat_obj_nds)), #fdc function in flacco
            domi = rename_domi(flacco::calculateFeatureSet(feat_obj_nds, set = "pca"))
          )
          
          # Save to CSV (header only for the first row)
          write.table(t(feature_vector), file = "s_ela_features.csv",
                      col.names = (write_count == 0), row.names = FALSE, append = TRUE, sep = ",")
          write_count <- write_count + 1
        })
        
        # Log progress
        cat(sprintf("Instance: %s | Repetition: %d | Time: %.2f sec\n",
                    instance_name, repetition_id, feature_time["elapsed"]))
      }
    }
  }
})

# Log total runtime
cat(sprintf("Total runtime: %.2f sec\n", total_time["elapsed"]))

