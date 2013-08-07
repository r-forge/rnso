isok <- function(x, bounds) {
  il <- min(x >= bounds[, 1])
  iu <- min(x <= bounds[, 2])
  isok <- min(il, iu)
  lisok <- (isok == 1)
  lisok
}