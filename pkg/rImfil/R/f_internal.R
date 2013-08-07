f_internal <- function(x, h, core_data){
  options                 <- core_data$options
  imfil_scale_aware       <- options$scale_aware
  imfil_noise_aware       <- options$noise_aware
  imfil_least_squares     <- options$least_squares
  imfil_simple_function   <- options$simple_function
  imfil_extra_argument    <- options$extra_argument
  exarg                   <- options$extra_arg_value
  fun_data 		   <- core_data$fun_data
  fn    		   <- fun_data$fn
  qbounds		   <- fun_data$qbounds
  dbounds		   <- fun_data$dbounds
  
  mx <- nrow(x)
  nx <- ncol(x)
  z <- x
  for (ix in 1:nx){
    z[, ix] <- dbounds*x[, ix] + qbounds[, 1]
  }
  func_type <- imfil_noise_aware + 10*imfil_scale_aware + 
		100*imfil_simple_function + 1000*imfil_extra_argument
  if (func_type == 0) {
    tmp <- fn(z)
    fx  <- tmp$fx
    iff <- tmp$iff
    icf <- tmp$icf
    tol <- 0
  } else     if (func_type == 10) {
    tmp <- fn(z, h)
    fx  <- tmp$fx
    iff <- tmp$iff
    icf <- tmp$icf
    tol <- 0
  } else   if (func_type == 1) {
    tmp <- fn(z)
    fx  <- tmp$fx
    iff <- tmp$iff
    icf <- tmp$icf
    tol <- tmp$tol
  } else   if (func_type == 11) {
    tmp <- fn(z, h)
    fx  <- tmp$fx
    iff <- tmp$iff
    icf <- tmp$icf
    tol <- tmp$tol
  } else   if (func_type == 100) {
    tmp <- fn(z)
    fx  <- tmp$fx
    mz  <- nrow(z)
    nz  <- ncol(z)
    iff <- matrix(0, nz, 1)
    icf <- matrix(nz, nz, 1)
    tol <- 0
  } else   if (func_type == 110) {
    tmp <- fn(z, h)
    fx  <- tmp$fx
    mz  <- nrow(z)
    nz  <- ncol(z)
    iff <- matrix(0, nz, 1)
    icf <- matrix(nz, nz, 1)
    tol <- 0
  } else   if (func_type == 101) {
    tmp <- fn(z, h)
    fx  <- tmp$fx
    mz  <- nrow(z)
    nz  <- ncol(z)
    iff <- matrix(0, nz, 1)
    icf <- matrix(nz, nz, 1)
    tol <- 0
  }  else   if (func_type == 111) {
    tmp <- fn(z, h)
    mz  <- nrow(z)
    nz  <- ncol(z)
    iff <- matrix(0, nz, 1)
    icf <- matrix(nz, nz, 1)
    tol <- 0
  }  else   if (func_type == 1000) {
    tmp <- fn(z, exarg)
    fx  <- tmp$fx
    iff <- tmp$iff
    icf <- tmp$icf
    tol <- 0
  } else   if (func_type == 1010) {
    tmp <- fn(z, h, exarg)
    fx  <- tmp$fx
    iff <- tmp$iff
    icf <- tmp$icf
    tol <- 0
  } else   if (func_type == 1001) {
    tmp <- fn(z, exarg)
    fx  <- tmp$fx
    iff <- tmp$iff
    icf <- tmp$icf
    tol <- tmp$tol
  } else   if (func_type == 1100) {
    tmp <- fn(z, exarg)
    fx  <- tmp$fx
    mz  <- nrow(z)
    nz  <- ncol(z)
    iff <- matrix(0, nz, 1)
    icf <- matrix(nz, nz, 1)
    tol <- 0
  } else { 
    stop("f_internal: func_type not found")
  }
  mf <- nrow(fx)
  nf <- ncol(fx)
  
  if (imfil_fscale <= 0){
    if(imfil_fscale == 0){
      imfil_fscale <- -1.2
    }
    if (imfil_least_squares == 1){
      val <- fx[, 1]*fx[, 1]/2
      scale_base <- val
    }else imfil_fscale <- abs(imfil_fscale)*scale_base
  }
  if (imfil_least_squares == 0){
    fx <- fx/imfil_fscale
  } else {
    fx <- fx/sqrt(imfil_fscale)
  }
  tol <- tol/imfil_fscale
  
  list(fx = fx, iff = iff, icf = icf, tol = tol)
  
}