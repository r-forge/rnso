imfil_inner_iteration <- function (fn, x, fx, sdiff, xc, gc, iteration_data, hess, fcount) {
  
  h <- iteration_data$h
  options <- iteration_data$options
  core_data <- iteration_data$core_data
  imfil_least_squares <- options$least_squares
  imfil_complete_history <- options$complete_history
  obounds <- iteration_data$obounds
  imfil_exec <- options$executive
  
  if (imfil_exec == 1) {
    imfil_exec_function <- options$executive_function
  }
  if (imfil_exec == 0) {
    if (imfil_least_squares == 0) {
      tmp <- imfil_qn_update(fn, x, fx, sdiff, xc, gc, iteration_data, hess)
      xp <- tmp$xp
      fval <- tmp$fvalp
      funs <- tmp$funs
      fct <- tmp$qfct
      iarm <- tmp$iarm
      diff_hist <- tmp$diff_hist
      nfail <- tmp$nfail
      hess <- tmp$hess
    } else {
      tmp <- imfil_gauss_newton(fn, x, fx, sdiff, xc, gc, iteration_data, hess)
      xp <- tmp$xp
      fval <- tmp$fvalp
      funs <- tmp$funs
      fct <- tmp$qfct
      iarm <- tmp$iarm
      diff_hist <- tmp$diff_hist
      nfail <- tmp$nfail
      hess <- tmp$hess
    }
  } else {
    tmp <- imfil_exec_function(fn, x, fx, sdiff, xc, gc, iteration_data, hess)
    xp <- tmp$xp
      fval <- tmp$fvalp
      funs <- tmp$funs
      fct <- tmp$qfct
      iarm <- tmp$iarm
      diff_hist <- tmp$diff_hist
      nfail <- tmp$nfail
      hess <- tmp$hess
  }
  fcount <- fcount + fct
  iteration_datap <- iteration_data
  if (imfil_complete_history > 0) {
    iteration_datap$complete_history <- many_point_hist_update(iteration_data$complete_history, diff_hist)
  }
  tmp <- reconcile_best_point(funs, xp, iteration_datap)
  iteration_datap <- tmp$new_data
  rflag <- tmp$rflag
  
  list(xp = xp, fval = fval, funs = funs, fcount = fcount, iarm = iarm,
  iteration_datap = iteration_datap, nfail = nfail, hess = hess)  
}