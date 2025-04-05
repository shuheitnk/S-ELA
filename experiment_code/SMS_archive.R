source("SMS-EMOA.R")
source("compute_hv.R")
biobj_info <- read.csv("bbob_biobj_ideal_nadir.csv", fileEncoding = "UTF-8")


# Define common options
common_option = list(
  optparse::make_option("--n.dim", type = "numeric", default = 2, help = "The dimensionality"),
  optparse::make_option("--n.obj", type = "numeric", default = 2, help = "The number of objectives"),
  optparse::make_option("--common_seed", type = "numeric", default = 2025L, help = "Common seed for reproducibility"),
  optparse::make_option("--n.repeat", type = "numeric", default = 31, help = "The number of repetitions"),
  optparse::make_option("--budget", type = "numeric", default = 20000L, help = "The maximum number of allowed function evaluations")
)

common_parser = optparse::OptionParser(option_list=common_option)
common_opt = optparse::parse_args(common_parser)

# Define SMS options
sms_option = list(
  optparse::make_option("--mu", type = "numeric", default = 100L),
  optparse::make_option("--mutator", type = "character", default = "mutPolynomial", help = ""),
  optparse::make_option("--mutPolynomial_p", type = "numeric", default = 0.2, help = "The probability of mutation"),
  optparse::make_option("--mutPolynomial_eta", type = "numeric", default = 20, help = "The distribution index of mutation"),
  optparse::make_option("--recombinator", type = "character", default = "recSBX", help = ""),
  optparse::make_option("--recSBX_p", type = "numeric", default = 0.9, help = "The probability of crossover"),
  optparse::make_option("--recSBX_eta", type = "numeric", default = 20, help = "The distribution index of crossover"),
  optparse::make_option("--log", type = "logical", default = TRUE, help = "Archive all solutions")
)

sms_parser = optparse::OptionParser(option_list=sms_option)
sms_opt = optparse::parse_args(sms_parser)

# Set seed and prepare unique random seeds for each repetition
set.seed(common_opt$common_seed)
random_seeds <- sample(1:100, common_opt$n.repeat, replace = FALSE)

# Loop over problem instances (iid)
for (iid in 1:5) {
  rep_id <- 0
  
  for (seed in random_seeds) {
    rep_id <- rep_id + 1
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
      
      # Run SMS-EMOA
      time_sms <- system.time({
        solutions <- get_sms_fitness(
          seed = seed,
          fitness.fun = fn,
          n.objectives = common_opt$n.obj,
          n.dim = common_opt$n.dim,
          minimize = TRUE,
          lower = fn_lower,
          upper = fn_upper,
          mu = sms_opt$mu,
          lambda = 1,
          ref.point = NULL,
          mutator = setup(mutPolynomial, p = sms_opt$mutPolynomial_p, eta = sms_opt$mutPolynomial_eta, lower = fn_lower, upper = fn_upper),
          recombinator = setup(recSBX, eta = sms_opt$recSBX_eta, p = sms_opt$recSBX_p, lower = fn_lower, upper = fn_upper),
          terminators = list(stopOnEvals(max.evals = common_opt$budget)),
          log.pop = sms_opt$log
        )
      })
      
      # Compute hypervolume
      nadir <- get_nadir(biobj_info, fid, iid)
      ideal_hv <- get_ideal_hv(biobj_info, fid, iid)
      non_dominated <- unique(solutions[, solutions[1, ] <= nadir[1] & solutions[2, ] <= nadir[2]])
      hv <- compute_hv(non_dominated, nadir)
      normalized_hv <- hv / ideal_hv
      
      # Save HV result
      hv_entry <- c(instance = fn_name, seed = seed, rep = rep_id, normalized_hv = normalized_hv)
      all_hv_data <- rbind(all_hv_data, hv_entry)
      
      # Print progress
      cat(sprintf("Instance: %s | Repetition: %d | Normalized_HV: %f | Elapsed time: %.2f sec\n", fn_name, rep_id, normalized_hv, time_sms["elapsed"]))
    }
  }
}
