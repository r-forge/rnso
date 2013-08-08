parallel_armijo <- function(fn, sdir, fold, xc, h, obounds, core_data) {
  options <- core_data$options
  imfil_least_squares <- options$least_squares
  imfil_limit_quasi_newton <- options$limit_quasi_newton
  imfil_verbose <- options$verbose
  beta <- options$armijo_reduction
  maxitarm <- options$maxitarm
  
  lambda <- 1
  n <- length(xc)
  fct <- 0
  aflag <- 1
  dd <- sdir
  fres <- matrix()
  frest <- fres
  diff_hist <- list(good_points = c(), good_values = c(), failed_points = c())
  if (imfil_limit_quasi_newton == 1) {
    smax <- 10*min(h, 1)
    if (sqrt(sum(d*d)) > smax) dd <- smax*dd/(sqrt(sum(dd)))
  }
  x <- xc
  fval <- fold
  
  number_steps <- maxitarm+1
  ddm <- matrix(0, n, number_steps)
  for(i in 1:number_steps) {
    ddm[, i] <- xc-lambda*dd
    lambda <- beta*lambda
    dd[, i] <- kk_proj(ddm[, i], obounds)
  }
  tmp <- fn(ddm, h, core_data)
  fta <- tmp$fv
  iflaga <- tmp$ifail
  ictra <- tmp$icount
  
  ibad <- (iflaga[, number_steps] == 1)
  if (sum(ibad) > 0) {
    diff_hist$failed_points <- ddm[, ibad]
  }
  igood <- (iflaga[, number_steps] == 0)
  if (sum(igood) > 0) {
    diff_hist$good_points <- ddm[, igood]
    if (imfil_least_squares == 1) {
      diff_hist$good_values <- fta[, igood]
    } else {
      diff_hist$good_values <- fta[, igood]
    }
  }
  if (imfil_least_squares == 1) {
    frest <- fta
  }
  
  fta <- f_to_vals(fta, imfil_least_squares)
  fct <- fct + sum(ictra)
  ilose <- (iflaga == 1)
  fta[ilose] <- fold + abs(fold)*0.0001 + 1e-12
  
  traditional <- 0
  if (traditional == 0) {
    ft <- min(fta)
    it <- which.min(fta)
    xt <- ddm[, it]
    iarm <- maxitarm
    if (ft < fval) {
      aflag <- 0
      fval <- ft
      if (imfil_least_squares == 1) fres <- frest[, it]
      x <- xt
      iarm <- min(it, maxitarm-1)
    } else {
      iarm <- -1
    }
    while(iarm < maxitarm && aflag == 1) {
      ft <- fta[iarm + 2]
      xt <- ddm[, iarm + 2]
      if (ft < fval && aflag == 1) {
	aflag <- 0
	if (imfil_least_squares == 1){
	fres <- frest[, iarm + 2]
      }
      fval <- ft
      x <- xt
    }
    iarm <- iarm + 1
  }
  if (iarm == maxitarm && aflag == 1 && imfil_verbose == 1) {
    print(cbind(iarm, h))
  }
  nfail <- aflag
  list(fct = fct, x = x,fval = fval, iarm = iarm, aflag = aflag,
      fres = fres, diff_hist = diff_hist, nfail = nfail)
  
}