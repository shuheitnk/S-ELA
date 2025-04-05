# Load necessary files and data
source("compute_hv.R")
biobj_info <- read.csv("bbob_biobj_ideal_nadir.csv", fileEncoding = "UTF-8")

# Define general options
common_option <- list(
  optparse::make_option("--n.dim", type = "numeric", default = 2, help = "Dimensionality of the problem"),
  optparse::make_option("--n.obj", type = "numeric", default = 2, help = "Number of objectives"),
  optparse::make_option("--common_seed", type = "numeric", default = 2025L, help = "Seed for reproducibility"),
  optparse::make_option("--n.repeat", type = "numeric", default = 31, help = "Number of repetitions"),
  optparse::make_option("--budget", type = "numeric", default = 20000L, help = "Maximum function evaluations")
)

common_parser <- optparse::OptionParser(option_list = common_option)
common_opt <- optparse::parse_args(common_parser)

# Define MOLE-specific options
mole_option <- list(
  optparse::make_option("--max_local_sets", type = "numeric", default = 20000L),
  optparse::make_option("--epsilon_gradient", type = "numeric", default = 1e-5),
  optparse::make_option("--descent_direction_min", type = "numeric", default = 1e-8),
  optparse::make_option("--descent_step_min", type = "numeric", default = 1e-5),
  optparse::make_option("--descent_step_max", type = "numeric", default = 0.01),
  optparse::make_option("--descent_scale_factor", type = "numeric", default = 2),
  optparse::make_option("--descent_armijo_factor", type = "numeric", default = 1e-4),
  optparse::make_option("--descent_history_size", type = "numeric", default = 100),
  optparse::make_option("--descent_max_iter", type = "numeric", default = 1000),
  optparse::make_option("--explore_step_min", type = "numeric", default = 1e-4),
  optparse::make_option("--explore_step_max", type = "numeric", default = 1e-1),
  optparse::make_option("--explore_angle_max", type = "numeric", default = 20),
  optparse::make_option("--explore_scale_factor", type = "numeric", default = 2),
  optparse::make_option("--refine_after_nstarts", type = "numeric", default = 10),
  optparse::make_option("--refine_hv_target", type = "numeric", default = 2e-5)
)

mole_parser <- optparse::OptionParser(option_list = mole_option)
mole_opt <- optparse::parse_args(mole_parser)

# Set seed and prepare unique random seeds for each repetition
set.seed(common_opt$common_seed)
random_seeds <- sample(1:100, common_opt$n.repeat, replace = FALSE)

# Loop over problem instances (iid)
for (iid in 1:5) {
  repetition_id <- 0
  
  for (seed in random_seeds) {
    repetition_id <- repetition_id + 1
    all_hv_data <- matrix(NA, nrow = 0, ncol = 4)
    
    # Define 55 BBOB bi-objective functions
    problems <- lapply(1:55, function(fid) {
      list(
        fn = smoof::makeBiObjBBOBFunction(dimensions = common_opt$n.dim, fid = fid, iid = iid),
        func_id = fid
      )
    })
    
    # Evaluate each problem
    for (problem in problems) {
      fn <- problem$fn
      fid <- problem$func_id
      fn_name <- smoof::getName(fn)
      fn_lower <- smoof::getLowerBoxConstraints(fn)
      fn_upper <- smoof::getUpperBoxConstraints(fn)
      
      # Generate starting points
      n_starts <- 100
      starting_points <- do.call(rbind, lapply(1:n_starts, function(x) {
        moleopt::runif_box(fn_lower, fn_upper)
      }))
      
      # Run MOLE optimization
      time_mole <- system.time({
        mole_trace <- moleopt::run_mole(
          fn, starting_points,
          max_local_sets        = mole_opt$max_local_sets,
          epsilon_gradient      = mole_opt$epsilon_gradient,
          descent_direction_min = mole_opt$descent_direction_min,
          descent_step_min      = mole_opt$descent_step_min,
          descent_step_max      = sqrt(sum((fn_upper - fn_lower)^2)) / 100,
          descent_scale_factor  = mole_opt$descent_scale_factor,
          descent_armijo_factor = mole_opt$descent_armijo_factor,
          descent_history_size  = mole_opt$descent_history_size,
          descent_max_iter      = mole_opt$descent_max_iter,
          explore_step_min      = mole_opt$explore_step_min,
          explore_step_max      = mole_opt$explore_step_max,
          explore_angle_max     = mole_opt$explore_angle_max,
          explore_scale_factor  = mole_opt$explore_scale_factor,
          refine_after_nstarts  = mole_opt$refine_after_nstarts,
          refine_hv_target      = mole_opt$refine_hv_target,
          max_budget            = common_opt$budget,
          logging               = "none"
        )
        
        solutions <- t(mole_trace$sets[[1]]$obj_space)
      })
      
      # Compute hypervolume
      nadir <- get_nadir(biobj_info, fid, iid)
      ideal_hv <- get_ideal_hv(biobj_info, fid, iid)
      non_dominated <- unique(solutions[, solutions[1, ] <= nadir[1] & solutions[2, ] <= nadir[2]])
      hv <- compute_hv(non_dominated, nadir)
      normalized_hv <- hv / ideal_hv
      
      # Save HV result
      hv_entry <- c(instance = fn_name, seed = seed, repetition = repetition_id, normalized_hv = normalized_hv)
      all_hv_data <- rbind(all_hv_data, hv_entry)
      
      # Print progress
      cat(sprintf("Instance: %s | Repetition: %d | Normalized_HV: %f | Elapsed time: %.2f sec\n", fn_name, repetition_id, normalized_hv, time_mole["elapsed"]))
    }
  }
}
