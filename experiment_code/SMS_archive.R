# Load libraries
library(optparse)
library(smoof)
library(ecr)
library(ecr3vis)
library(plyr)

source("SMS-EMOA.R")
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

# Define SMS options
sms_option = list(
  make_option("--mu", type = "numeric", default = 100L),
  make_option("--mutator", type = "character", default = "mutPolynomial", help = ""),
  make_option("--mutPolynomial_p", type = "numeric", default = 0.2, help = "The probability of mutation"),
  make_option("--mutPolynomial_eta", type = "numeric", default = 20, help = "The distribution index of mutation"),
  make_option("--recombinator", type = "character", default = "recSBX", help = ""),
  make_option("--recSBX_p", type = "numeric", default = 0.9, help = "The probability of crossover"),
  make_option("--recSBX_eta", type = "numeric", default = 20, help = "The distribution index of crossover"),
  make_option("--log", type = "logical", default = TRUE, help = "Archive all solutions")
)

sms_parser = OptionParser(option_list=sms_option)
sms_opt = parse_args(sms_parser)

set.seed(common_opt$common_seed)
random_numbers_unique <- sample(1:100, common_opt$n.repeat, replace = FALSE)
random_numbers_unique <- random_numbers_unique[12]


# Loop through random seeds
for (iii in 10:10) {
  roop <- 11
  
  for (roop_seed in random_numbers_unique) {
    # Pre-allocate storage for Archive_Data and all_hv_data
    Archive_Data <- matrix(NA, nrow = 0, ncol = common_opt$budget + 4)
    all_hv_data <- matrix(NA, nrow = 0, ncol = 4)
    roop <- roop + 1
    
    
    # Define optimization problems
    problems <- lapply(43:55, function(fid) {
      list(fn = makeBiObjBBOBFunction(dimensions = 2, fid = fid, iid = iii), func_id = fid)
    })
    
    # Loop through each problem
    for (problem in problems) {
      obj.fn <- problem$fn
      func_id <- problem$func_id
      
      instance_name <- smoof::getName(obj.fn)
      fn.lower <- smoof::getLowerBoxConstraints(obj.fn)
      fn.upper <- smoof::getUpperBoxConstraints(obj.fn)
      
      time_sms <- system.time({
        # Run SMS-EMOA optimization
        sms_solution <- get_sms_fitness(
          seed = roop_seed,
          fitness.fun = obj.fn,
          n.objectives = common_opt$n.obj,
          n.dim = common_opt$n.dim,
          minimize = TRUE,
          lower = fn.lower,
          upper = fn.upper,
          mu = sms_opt$mu,
          lambda = 1,
          ref.point = NULL,
          mutator = setup(mutPolynomial, p = sms_opt$mutPolynomial_p, eta = sms_opt$mutPolynomial_eta, lower = fn.lower, upper = fn.upper),
          recombinator = setup(recSBX, eta = sms_opt$recSBX_eta, p = sms_opt$recSBX_p, lower = fn.lower, upper = fn.upper),
          terminators = list(stopOnEvals(max.evals = common_opt$budget)),
          log.pop = sms_opt$log
        )
      })
      
      nadir_point <- get_nadir(biobj_info, func_id, iii)
      better_sol <- unique(sms_solution[, sms_solution[1, ] <= nadir_point[1] & sms_solution[2, ] <= nadir_point[2]])
      hv <- compute_hv(better_sol, nadir_point)
      Ideal_hv <- get_ideal_hv(biobj_info, func_id, iii)
      
      normalized_hv <- hv / Ideal_hv
      hv_data <- c(instance = instance_name, seed = roop_seed, roop = roop, normalized_hv = normalized_hv)
      
      # Append HV data directly to the matrix
      all_hv_data <- rbind(all_hv_data, hv_data)
      
      name_list <- cbind(
        matrix(c(instance_name, instance_name), nrow = 2, ncol = 1),
        matrix(c("y1", "y2"), nrow = 2, ncol = 1),
        matrix(c(roop_seed, roop_seed), nrow = 2, ncol = 1),
        matrix(c(roop, roop), nrow = 2, ncol = 1),
        sms_solution
      )
      
      # Append the solution data to Archive_Data
      Archive_Data <- rbind(Archive_Data, name_list)
      
      print(paste("instance name ", instance_name, "roop: ", roop))
      print(paste("process time: ", time_sms["elapsed"], "sec"))
    }
    
    # Write HV data and solution data to file
    file_name_hv <- paste0("sms_hv_instance", iii, ".csv")
    write.table(all_hv_data, file = file_name_hv, col.names = FALSE, row.names = FALSE, append = TRUE, sep = ",")
    
    
    
    file_name <- paste0("sms_archive_instance", iii, ".csv")
    write.table(t(Archive_Data), file = file_name, col.names = FALSE, row.names = FALSE, append = TRUE, sep = ",")
    
  }
  
}
