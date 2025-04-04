compute_hv <- function (solution_set, reference_point){
  
  
  if (length(solution_set)<3){
    if (length(solution_set)==0){
      HV <- 0
    }else{
      HV <- (reference_point[1]-solution_set[1])*(reference_point[2]-solution_set[2])
    }
    
    
  }else{
    pareto.matrix <- solution_set[, ecr::nondominated(t(solution_set))]
    
    if(is.matrix(pareto.matrix)){
      HV <- ecr3vis::hv(pareto.matrix, reference_point)
    }else {
      HV <- 0
    }
    
  }
 
  return(HV)
}

get_ideal_hv <- function(biobj_info,f, i){
  
  instance_data <- subset(biobj_info, fid == f & iid == i)
  Ideal_hv <- instance_data$ideal_hv
  
  return(Ideal_hv)
}