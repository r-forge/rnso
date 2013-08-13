imfil_gauss_newton <- function (fn, x, fun, jac, xc,  gc, iteration_data, hessold) {
  obounds <- iteration_data$obounds
  options <- iteration_data$options
  core_data <- iteration_data$core_data
  h <- iteration_data$h
  
  funs <- fun
  fval <- sum(fun*fun)/2
  sgrad <- sum(jac*fun)
  
  hess <- sum(jac*jac)
  imfil_verbose <- options$verbose
  n <- length(x)
  
  epsb <- 1e-6
  alist <- max((x > obounds[, 1] + epsb), (x < obounds[, 2] - epsb))
  pr <- diag(alist)
  
  rjac <- jac%*%pr
  tmp <- qr(rjac)
  rq <- qr.Q(tmp)
  rr <- qr.R(tmp)
  sdir1 <- (diag(n) - pr) %*% sgrad
  sdir2 <- sum(rq*fun)
  sdir2 <- solve(rr, sdir2)
  sdir <- sdir1 + sdir2
  
  tmp <- armijo_explore(fn, sdir, fval, x, h, core_data, obounds)
  qfct <- tmp$fct
  xp <- tmp$x
  fvalp <- tmp$fval
  iarm <- tmp$iarm
  fres <- tmp$fres
  diff_hist <- tmp$diff_hist
  nfail <- tmp$nfail
  
  if (length(fres) > 0) {
    funs <- fres
  }
  list(xp = xp, fvalp = fvalp, funs = funs, qfct = qfct, iarm = iarm,
  diff_hist = diff_hist, nfail = nfail, hess = hess)
}
