jac_or_grad <- function(sgrad, jac, options) {
  imfil_least_squares <- options$least_squares
  if (imfil_least_squares == 1) {
    sdiff <- jac
  } else if (imfil_least_squares == 0) {
    sdiff <- sgrad
  } else {
    stop("jac_or_grad: wrong least_squares")
  }
  sdiff
}