

get_nadir <- function(biobj_info,f, i) {
  
  
  instance_data <- subset(biobj_info, fid == f & iid == i)
  nadir_point <- c(instance_data$nadir1,instance_data$nadir2)
  
  return(nadir_point)
  
}