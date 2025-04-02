tchebycheff <- function(Y, w){
  z <- apply(Y, 2, min)
  G <- sapply(1:nrow(Y), function(j) {
    max(w * (Y[j, ] - z))
  })
  return(G)
}
