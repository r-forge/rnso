random_augment <- function(vin, k){
  vin <- as.matrix(vin)
  n  <- nrow(vin)
  nn <- ncol(vin)
  if (k == 0) {
    vout <- vin
    return(vout)
  }
  rv <- matrix(runif(n*k), n, k)
  for (ir in 1:k){
    rv[, ir] <- rv[, ir]/sqrt(sum(rv[, ir]*rv[, ir]))
  }
  vout <- cbind(vin, rv)
  vout
}