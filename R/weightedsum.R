weightedsum <- function(Y, w){
  G <- sapply(1:nrow(Y), function(i) {
    w[i] * Y[i, j]
  })
  return(G)
}
