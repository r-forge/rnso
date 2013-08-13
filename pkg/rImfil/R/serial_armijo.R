serial_armijo <- function (fn, sdir, fold, xc, h, obounds, core_data) {
  options <- core_data$options
  imfil_least_squares <- options$least_squares
  imfil_limit_quasi_newton <- options$limit_quasi_newton
  imfil_verbose <- options$verbose
  beta <- options$armijo_reduction
  maxitarm <- options$maxitarm
  
  lambda <- 1
  n <- length(xc)
  iarm <- -1
  fct <- 0
  aflag <- 1
  dd <- sdir
  fres <- matrix()
  frest <- fres
  diff_hist <- list(good_points = matrix(), good_values = matrix(),
		failed_points = matrix())
  if (imfil_limit_quasi_newton == 1) {
    smax <- 10*min(h, 1) 
    if (sqrt(sum(dd)) > smax) {
      dd <- smax*dd/sqrt(sum(dd))
    }
  }
  x <- xc
  fval <- fold
  while (iarm < maxitarm && aflag == 1) {
    d <- -lambda*dd
    xt <- x + d
    xt <- kk_proj(xt, obounds)
    tmp <- fn(xt, h, core_data)
    ft <- tmp$fv
    ifl <- tmp$ifail
    ict <- tmp$icount
    
    fct <- fct + ict
    diff_hist <- single_point_hist_update(diff_hist, xt, ft, ifl)
    if (imfil_least_squares == 1) {
      frest <- ft
      ft <- sum(frest*frest)/2
    }
    if (ifl == 1) {
      ft <- fold + abs(fold)*0.001 + 1e-12
    }
    if (ft < fval && aflag == 1) {
      aflag <- 0 
      fval <- ft
      fres <- frest
      x <- xt
    }
    if (aflag == 1) {
      lambda <- beta*lambda
    }
    iarm <- iarm +1
  }
  
  if (iarm == maxitarm && aflag == 1 && imfil_verbose == 1) {
    cat ("line search failure", cbind(iarm, h))
  }
  nfail <- aflag
  list(fct = fct, x = x, fval = fval, iarm = iarm, aflag = aflag,
	fres = fres, diff_hist = diff_hist, nfail = nfail)  
}