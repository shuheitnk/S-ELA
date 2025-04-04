# Load libraries
library(optparse)
library(smoof)
library(ecr)
library(ecr3vis)
library(plyr)
library(moleopt)
library(purrr)

source("MakeBiObjBBOB.R")
source("compute_hv.R")
source("get_nadir.R")
biobj_info <- read.csv("bbob_biobj_ideal_nadir.csv", fileEncoding = "UTF-8")

# Define common options
common_option = list(
  make_option("--n.dim", type = "numeric", default = 2, help = "The dimensionality"),
  make_option("--n.obj", type = "numeric", default = 2, help = "The number of objectives"),
  make_option("--common_seed", type = "numeric", default = 2025L, help = "Common seed for reproducibility"),
  make_option("--n.repeat", type = "numeric", default = 31, help = "The number of repetitions"),
  make_option("--budget", type = "numeric", default = 20000L, help = "The maximum number of allowed function evaluations")
)

common_parser = OptionParser(option_list=common_option)
common_opt = parse_args(common_parser)

mole_list = list(
  make_option("--max_local_sets", type = "numeric", default = 20000L),
  make_option("--epsilon_gradient", type = "numeric", default = 1e-5, help = ""),
  make_option("--descent_direction_min", type = "numeric", default = 1e-8, help = ""),
  make_option("--descent_step_min", type = "numeric", default = 1e-5, help = ""),
  make_option("--descent_step_max", type = "numeric", default = 0.01, help = ""),
  make_option("--descent_scale_factor", type = "numeric", default = 2, help = ""),
  make_option("--descent_armijo_factor", type = "numeric", default = 1e-4, help = "[recCrossover, recIntermediate, recOX, recPMX, recSBX]"),
  make_option("--descent_history_size", type = "numeric", default = 100, help = ""),
  make_option("--descent_max_iter", type = "numeric", default = 1000, help = ""),
  make_option("--explore_step_min", type = "numeric", default = 1e-4, help = ""),
  make_option("--explore_step_max", type = "numeric", default = 1e-1, help = ""),
  make_option("--explore_angle_max", type = "numeric", default = 20, help = ""),
  make_option("--explore_scale_factor", type = "numeric", default = 2, help = ""),
  make_option("--refine_after_nstarts", type = "numeric", default = 10, help = ""),
  make_option("--refine_hv_target", type = "numeric", default = 2e-5, help = "")
)

mole_parser = OptionParser(option_list=mole_list)
mole_opt = parse_args(mole_parser)


set.seed(common_opt$common_seed)
random_numbers_unique <- sample(1:100, common_opt$n.repeat, replace = FALSE)
random_numbers_unique <- random_numbers_unique[1:common_opt$n.repeat]


# Loop through random seeds
for (iid in 1:5) {
  roop <- 0
  
  for (roop_seed in random_numbers_unique) {
    # Pre-allocate storage for Archive_Data and all_hv_data
    
    all_hv_data <- matrix(NA, nrow = 0, ncol = 4)
    roop <- roop + 1
    
    
    # Define optimization problems
    problems <- lapply(1:55, function(fid) {
      list(fn = makeBiObjBBOBFunction(dimensions = 2, fid = fid, iid = iid), func_id = fid)
    })
    
    # Loop through each problem
    for (problem in problems) {
      obj.fn <- problem$fn
      func_id <- problem$func_id
      
      instance_name <- smoof::getName(obj.fn)
      fn.lower <- smoof::getLowerBoxConstraints(obj.fn)
      fn.upper <- smoof::getUpperBoxConstraints(obj.fn)
      
      
      nstarts <- 100
      starting_points <- lapply(1:nstarts, function(x) runif_box(fn.lower, fn.upper))
      starting_points <- do.call(rbind, starting_points)
      
      time_mole <- system.time({
        # Run mole optimization
        
        mole_trace <- run_mole(obj.fn, starting_points,
                               max_local_sets = 1000,
                               epsilon_gradient = 1e-8,
                               descent_direction_min = 1e-6,
                               descent_step_min = 1e-6,
                               descent_step_max = sqrt(sum((fn.upper - fn.lower) ** 2)) / 100,
                               descent_scale_factor = 2,
                               descent_armijo_factor = 1e-4,
                               descent_history_size = 100,
                               descent_max_iter = 1000,
                               explore_step_min = 1e-4,
                               explore_step_max = 1e-2,
                               # explore_step_max = sqrt(sum((upper - lower) ** 2)) / 100,
                               explore_angle_max = 20,
                               explore_scale_factor = 2,
                               refine_after_nstarts = 100,
                               refine_hv_target = 2e-5,
                               # custom_descent_fn = create_lbfgsb_descent(f, lower, upper),
                               # lower = rep(-100, length(lower)),
                               # upper = rep(100, length(upper)),
                               max_budget = common_opt$budget,
                               logging = "none"
        )
        
        mole_solution <- t(mole_trace$sets[[1]]$obj_space)
        
      })
      
      nadir_point <- get_nadir(biobj_info, func_id, iii)
      better_sol <- unique(mole_solution[, mole_solution[1, ] <= nadir_point[1] & mole_solution[2, ] <= nadir_point[2]])
      hv <- compute_hv(better_sol, nadir_point)
      
      Ideal_hv <- get_ideal_hv(biobj_info, func_id, iii)
      
      normalized_hv <- hv / Ideal_hv
      hv_data <- c(instance = instance_name, seed = roop_seed, roop = roop, normalized_hv = normalized_hv)
      
      
      # Append HV data directly to the matrix
      all_hv_data <- rbind(all_hv_data, hv_data)
      
      print(paste("instance name ", instance_name, "roop: ", roop))
      print(paste("process time: ", time_mole["elapsed"], "sec"))
    }
    
    # Write HV data and solution data to file
    #file_name_hv <- paste0("mole_hv_instance", iid, ".csv")
    #write.table(all_hv_data, file = file_name_hv, col.names = FALSE, row.names = FALSE, append = TRUE, sep = ",")
    
  }
  
}
