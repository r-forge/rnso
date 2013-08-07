scan_history <- function(complete_history, xp, fp, dx) {
  mp <- nrow(xp)
  np <- ncol(xp)
  imold <- matrix(0, 1, np)
  oldpoints <- c()
  oldvalues <- c()
  oldflags  <- c()
  for (i in 1:np) {
    tmp <- scal_complete_history(complete_history, xp[, i])
    fpt <- tmp$fhist
    ift <- tmp$iflaghist
    if (ift > -1) {
      oldpoints <- cbind(oldpoints, xp[, i])
      oldvalues <- cbind(oldvalues, fpt)
      oldflags  <- cbind(oldflags, ift)
      imold[i]  <- 1
    }
  }
  oldindex <- (imold == 1)
  list(oldindex = oldindex, oldpoints = oldpoints, oldvalues = oldvalues,
	oldflags = oldflags)
}