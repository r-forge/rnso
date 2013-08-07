f_to_vals <- function(funs, least_squares){
  m <- nrow(funs)
  n <- ncol(funs)
  if (m == 0 && n == 0){
    fvals <- c()
  }
  if(least_squares == 1){
    for(i in 1:n){
      fvals[i] <- sum(funs[, i]*funs[, i])/2
    }
  } else {
    fvals = funs
  }
  fvals
}