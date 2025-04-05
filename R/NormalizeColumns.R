NormalizeColumns <- function(Y) {

  if (nrow(Y) == 0 || ncol(Y) == 0) {
    stop("Input matrix is empty.")
  }

  Y_min <- apply(Y, 2, min)
  Y_max <- apply(Y, 2, max)
  Y_normalized <- sweep(sweep(Y, 2, Y_min, FUN = "-"), 2, Y_max - Y_min, FUN = "/")

  return(Y_normalized)
}
