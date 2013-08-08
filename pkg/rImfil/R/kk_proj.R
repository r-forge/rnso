kk_proj <- function(x, bounds) {
  ndim <- length(x)
  px <- matrix(0, ndim, 1)
  px <- min(bounds[, 2], x)
  px <- max(bounds[, 1], px)
  px
}