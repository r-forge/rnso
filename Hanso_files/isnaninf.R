is.nainf <- function(v) {
  any(is.na(v)) || any(is.infinite(v)) || any(is.nan(v))
}