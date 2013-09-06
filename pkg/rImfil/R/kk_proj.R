kk_proj <- function(x, bounds) {
  ndim <- length(x)
  px <- matrix(0, ndim, 1)
  px <- apply(cbind(bounds[, 2],x), 1,min)
  px <- apply(cbind(bounds[, 1],x), 1,min)
  px
}