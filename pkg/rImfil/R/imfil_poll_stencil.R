imfil_poll_stencil <- function(x, fn, dx, dc, bounds, core_data, h, complete_history) {
  options       <- core_data$options
  parallel      <- options$parallel
  least_squares <- options$least_squares
  n <- nrow(dx)
  vsize <- ncol(dx)
  
  best_point <- x
  if (least_squares == 1){
    best_value <- sum(fc*fc)/2
  } else {
    best_value <- fc
  } 
  best_value_f <- fc
  
  iflag <- matrix(0, vsize, 1)
  m <- length(fc)
  fp <- matrix(0, m, vsize)
  failed_points <- c()
  good_points   <- c()
  good_points   <- c()
  good_values   <- c()
  sgood <- 0
  fval <- matrix(0, vsize, 1)
  icount <- 0
  pold <- 0
  dx1 <- c()
  xp1 <- c()
  stencil_type <- options$stencil
  if (stencil_type == 1) {
    for (i in 1:vsize) {
      xp[, 1] <- x + dx[, i]
      if (isok(xp[, i], bounds) == 0) {
	dx[, i] <- -dx[, i]
      }
    }
  }
  for (i in 1:vsize) {
    xp[, i] <- x + dx[, i]
    if (isok(xp[, i], bounds)){
      pold <- pold +1
      dx1 <- cbind(dx1, dx[, i])
      xp1 <- cbind(xp1, xp[, i])
    }
  }
  fp1 <- fp[, 1:pold]
  xp <- xp1
  dx <- dx1
  fp <- fp1
  tmp <- scan_history(complete_history, xp, fp, dx)
  oldindex   <- tmp$oldindex
  oldpoints  <- tmp$oldpoints
  oldvalues  <- tmp$oldvalues
  oldflags   <- tmp$oldflags
  newindex <- !oldindex
  xp <- xp1[, newindex]
  iflago <- matrix(0, 1, vsize)
  if (sum(oldindex) > 0) {
    fp[, oldindex] <- oldvalues
    iflago[oldindex] <- oldflags
  }
  pnew <- sum(newindex)
  fp1 <- c()
  iflag <- c()
  if (parallel == 0) {
    for (i in 1:pnew) {
      tmp <- fn(xp, h, core_data)
      fpx    <- tmp$fv
      iflag  <- tmp$ifail
      icount <- tmp$ict
      fp1 <- cbind(fp1, fpx)
      iflag <- cbind(iflag, iflagx)
      icount <- icount + ict
    }
    } else {
      if (pnew > 0) {
	tmp <- fn(xp, h, core_data)
	fp1    <- tmp$fv
	iflag  <- tmp$ifail
	ictrp  <- tmp$ict
	icount <- icount + sum(ictrp) 
      }
    }  
  if (pnew > 0) {
    fp[, newindex] <- fp1
    iflago[newindex] <- iflag
  }
  fp1 <- fp
  iflag <- iflago
  ibad <- (iflag[1:pold] == 1)
  if (sum(ibad) > 0) {
    failed_points <- xp1[, ibad]
  }
  igood <- (iflag[1:pold] == 0)
  sgood <- sum(igood)
  good_dx <- c()
  good_df <- c()
  if (sgood > 0){
    good_points <- xp1[, igood]
    if (least_squares == 1) {
      good_fp <- fp1[, igood]
    } else {
      good_fp <- fp1[igood]
    }
    good_dx <- dx1[, igood]
    for (ig in 1:sgood){
      if (least_squares == 1){
	good_df[, ig] <- good_fp[, ig] - fc
      } else {
	good_df[ig] <- good_fp[ig] - fc
      }
    }
    good_values <- good_fp
  }
  list(best_value = best_value, best_value_f = best_value_f, best_point = best_point,
      icount = icount, sgood = sgood, good_points = good_points, good_values = good_values,
      good_dx = good_dx, good_df = good_df, failed_points = failed_points)
}