armijo_explore <- function (fn, sdir, fold, xc, h, core_data, obounds) {
  options <- core_data$options
  imfil_parallel <- options$parallel
  imfil_least_squares <- options$least_squares
  imfil_verbose <- options$verbose
  beta <- imfil$armijo_reduction
  maxitarm <- imfil$maxitarm
  
  if (imfil_parallel == 0) {
    tmp <- serial_armijo(fn, sdir, fold, xc, h, obounds, core_data)
    fct <- tmp$fct
    x <- tmp$x
    fval <- tmp$fval
    iarm <- tmp$iarm
    aflag <- tmp$aflag
    fres <- tmp$fres
    diff_hist <- tmp$diff_hist
    nfail <- tmp$nfail
  } else {
    tmp <- parallel_armijo(fn, sdir, fold, xc, h, obounds, core_data)
    fct <- tmp$fct
    x <- tmp$x
    fval <- tmp$fval
    iarm <- tmp$iarm
    aflag <- tmp$aflag
    fres <- tmp$fres
    diff_hist <- tmp$diff_hist
    nfail <- tmp$nfail
  }
  if (iarm == maxitarm && aflag == 1 && imfil_verbose == 1) {
    cat ("line search failure", cbind(iarm, h))
  }
  list(fct = fct, x = x, fval = fval, iarm = iarm, aflag = aflag,
	fres = fres, diff_hist = diff_hist, nfail = nfail)  
}