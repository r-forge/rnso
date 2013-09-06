imfil_qn_update <- function (fn, x, fval, sgrad, xc, gc, iteration_data,
  hessold) {
  
  core_data <- iteration_data$core_data
  obounds <- iteration_data$obounds
  h <- iteration_data$h
  
  itc <- iteration_data$itc
  options <- core_data$options
  imfil_verbose <- options$verbose
  quasi <- options$quasi
  
  nx <- length(x)
  hess <- hessold
  
  if (itc > 1) {
    epsb <- 1e-6
    alist <- (x > obounds[, 1] + epsb) && (x < obounds[, 2] - epsb)
    if (quasi == 1) {
      hess <- bfupdate(x, xc, sgrad, gc, hessold, alist)
    } else if (hess == 1) {
      hess <- sr1up(x, xc, sgrad, gc, hessold, alist)
    } else {
      hess <- diag(nx)
    }
  }
  sdir <- solve(hess, sgrad)
  tmp <- armijo_explore(fn, sdir, fval, x, h, core_data, obounds)
  qfct <- tmp$fct
  xp <- tmp$x
  fvalp <- tmp$fval
  iarm <- tmp$iarm
  fres <- tmp$fres
  diff_hist <- tmp$diff_hist
  nfail <- tmp$nfail
  
  funs <- fvalp
  
  list(xp = xp, fvalp = fvalp, funs = funs, qfct = qfct, iarm = iarm,
  diff_hist = diff_hist, nfail = nfail, hess = hess)
}