scal_complete_history <- function(complete_history, x){
  iflaghist <- -1
  losers <- complete_history$failed_points
  losers = as.matrix(losers)
  ml <- nrow(losers)
  nl <- ncol(losers)
  winners <- complete_history$good_points
  mw <- nrow(winners)
  nw <- ncol(winners)
  fhist <- matrix(NA, mw, 1)
  iquit <- 0
  if (!is.na(winners)){
    for (i in 1:nw){
      d <- norm(as.matrix(x-winners[, i]), "I")
      if (d < 1e-12 ){
        iquit <- 1
        fhist <- complete_history$good_values[, i]
        iflaghist <- 0
        break
      }
    }
  }
  browser()
  if (iquit == 0 && !is.null(losers) && !is.na(nl)) {
    for (i in 1:nl) {
      d <- norm(as.matrix(x-losers[, i]), 'I')
      if (d < 1e-12) {
	fhist <- NA
	iquit <- 1
	iflaghist <- 1
	break
      }
    }
  }
  list(fhist = fhist,iflaghist = iflaghist)
}