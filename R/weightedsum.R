weightedsum <- function(Y, w){
  G <- sapply(1:nrow(Y), function(j) {
    sum(w * Y[j, ])
  })
  return(G)
}