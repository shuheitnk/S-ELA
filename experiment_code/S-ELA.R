
# Source necessary files

source("calculate_new_feature.R")
source("makeObject.R")


# Initialize parameters
d = 2
set.seed(100)
uni <- sample(1:1000, 1, replace = FALSE)
ite = 1:1

nnn <- 0

# Measure the total time for execution
time_all <- system.time({
  
  ins = 0
  for (iid in ite) {
    
    ins <- ins + 1
    problems <- lapply(1:1, function(fid) {
      list(fn = smoof::makeBiObjBBOBFunction(dimensions = 2, fid = fid, iid = iid))
    })
    
    # Loop through each problem
    for (problem in problems) {
      roop = 0
      
      # Loop through the unique samples
      for (u in uni) {
        roop <- roop + 1
        
        time_feat <- system.time({
          
          # Set seed for reproducibility
          set.seed(u)
          
          # Get the objective function and instance name
          obj.fn <- problem$fn
          instance_name <- smoof::getName(obj.fn)
          
          # Generate samples and compute objective values
          Obj = makeobject(fn = obj.fn, n.sample = 400)
          X = Obj$sample
          Y = Obj$y
          obj_scale = Obj$obj_scale
          # Rank the objective values
          ranks <- ecr::doNondominatedSorting(t(Y))$ranks
          ranks <- as.numeric(ranks)
          
          r <- length(ranks)
          rank1_ratio <- sum(ranks == 1) / r
          rank_count <- table(ranks)
          
          # Preinitialize vectors for calculations
          weight <- MOEADr::decomposition_sld(list(name = "sld", H = 2, .nobj = 2))[, 1]
          n <- length(weight)
          
          feat_object_dom <- flacco::createFeatureObject(X = X, y = ranks)
          
          
          # Initialize result lists
          results_meta <- vector("list", n)
          results_distr <- vector("list", n)
          results_disp <- vector("list", n)
          results_nbc <- vector("list", n)
          results_ic <- vector("list", n)
          results_fdc <- vector("list", n)
          results_pca <- vector("list", n)
          
          # Loop through each weight and calculate features
          for (i in seq_along(weight)) {
            w <- weight[i]
            
            # Calculate new objective values based on the weight
            Fn <- sapply(1:nrow(Y), function(j) {
              w * Y[j, 1] + (1 - w) * Y[j, 2]
            })
            
            feat_object <- flacco::createFeatureObject(X = X, y = Fn)
            
            # Calculate various features
            results_meta[[i]] <- flacco::calculateFeatureSet(feat_object, set = "ela_meta")
            results_distr[[i]] <- flacco::calculateFeatureSet(feat_object, set = "ela_distr")
            results_disp[[i]] <- flacco::calculateFeatureSet(feat_object, set = "disp")
            results_nbc[[i]] <- flacco::calculateFeatureSet(feat_object, set = "nbc")
            results_ic[[i]] <- flacco::calculateFeatureSet(feat_object, set = "ic")
            results_fdc[[i]] <- flacco:::calculateFitnessDistanceFeatures (feat_object)
            results_pca[[i]] <- flacco::calculateFeatureSet(feat_object, set = "pca")
            
            
          }
          
          # Create statistics list
          stats_list <- c(
            instance = instance_name,
            obj_scale = minmax(obj_scale),
            rank1_ratio = rank1_ratio,
            
            # decomposition based features
            deco_meta = unlist(results_arrange(results_meta),, use.names = TRUE),
            deco_distr = unlist(results_arrange(results_distr)),
            deco_disp = unlist(results_arrange(results_disp)),
            deco_nbc = unlist(results_arrange(results_nbc)),
            deco_ic = unlist(results_arrange(results_ic)),
            deco_fdc = unlist(results_arrange(results_fdc)),
            deco_pca = unlist(results_arrange(results_pca)),
            
            # non-dominated sorting based features
            domi_meta = unlist(flacco::calculateFeatureSet(feat_object_dom, set = "ela_meta")),
            domi_distr = unlist(flacco::calculateDistributionFeatures(ranks)),
            domi_disp = unlist(flacco::calculateFeatureSet(feat_object_dom, set = "disp")),
            domi_nbc = unlist(flacco::calculateFeatureSet(feat_object_dom, set = "nbc")),
            domi_ic = unlist(flacco::calculateFeatureSet(feat_object_dom, set = "ic")),
            domi_fdc = unlist(FitnessDistanceFeatures(X = X, Y = ranks), use.names = TRUE),
            domi_pca = unlist(flacco::calculateFeatureSet(feat_object_dom, set = "pca"))
            
            
            
            
          )
          
          
          # Write results to CSV
          if (nnn == 0) {
            #write.table(t(stats_list), file = "test_fest_400.csv", col.names = TRUE, row.names = FALSE, append = TRUE, sep = ",")
          } else {
            #write.table(t(stats_list), file = "test_fest_400.csv", col.names = FALSE, row.names = FALSE, append = TRUE, sep = ",")
          }
          
          # Increment counter
          nnn <- nnn + 1
          
        })
        
        # Print progress
        print(paste("Instance name: ", instance_name, " - Roop: ", roop))
        print(paste("Execution time: ", time_feat["elapsed"], "sec"))
      }
      
    }
    
  }
  
})

# Print total execution time
print(paste("Total execution time: ", time_all["elapsed"], "sec"))
