stencil_diff <- function(x, fn, dx, fc, iteration_data, complete_history) {
  core_data <- iteration_data$core_data
  h <- iteration_data$h
  bounds <- iteration_data$bounds
  options <- core_data$options
  parallel <- options$parallel
  least_squares <- options$least_squares
  #browser()
  tmp <- imfil_poll_stencil(x, fn, dx, fc, bounds, core_data, h, complete_history)
  best_value <- tmp$best_value
  best_value_f <- tmp$best_value_f
  best_point <- tmp$best_point
  icount <- tmp$icount
  sgood <- tmp$sgood
  good_points <- tmp$good_points
  good_values <- tmp$good_values
  good_dx <- tmp$good_dx
  good_df <- tmp$good_df
  failed_points <- tmp$failed_points
  
  if (least_squares == 1) {
    fval <- sum(fc*fc)/2
  } else fval <- fc
  #browser()
  tmp <- collect_stencil_data(good_points,good_values,failed_points, 
                x,fval,fc,fc,options)
  sflag        <- tmp$sflag
  best_value   <- tmp$best_value
  best_value_f <- tmp$best_value_f
  best_point   <- tmp$best_point
  svar         <- tmp$svar
  diff_hist    <- tmp$diff_hist
  
  if (sgood > 0) {
    pdx <- ginv(good_dx)
    fprime <- good_df%*%pdx
    if (least_squares == 0) {
      grad <- t(fprime)
      jac <- matrix()
    } else {
      jac <- fprime
      grad <- t(fprime)%*%fc
    }
  } else {
    grad <- NA*x
    jac <- matrix()
  }
  
  list(grad = grad, best_value = best_value, best_value_f = best_value_f,
  best_point = best_point, icount = icount, sflag = sflag, svar = svar,
  diff_hist = diff_hist, jac = jac)
}