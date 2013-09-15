collect_stencil_data <- function (good_points, good_values, failed_points,
	best_point_old, best_value_old, best_value_f_old, fc, options) {
 
  browser()
  least_squares <- options$least_squares
  good_scalers <- f_to_vals(good_values, least_squares)
  xbest_value <- min(good_scalers)
  ibest <- min(which.min(good_scalers),1) # fix, sort of
  good_points <- as.matrix(good_points)
  xbest_point <- good_points[, ibest]
  xbest_value_f <- good_values[ibest] # warning! check
  best_value <- best_value_old
  best_value_f <- best_value_f_old
  best_point <- best_point_old
  browser()
  if(!is.na(xbest_value)){ #a quick fix, might slow down convergence.
    if (xbest_value < best_value_old) {
      best_value_old <- xbest_value
      best_value_old <- xbest_value_f
      best_point <- xbest_point
  }} else {
    best_value <- best_value_old
    best_value_f <- best_value_f_old
    best_point <- best_point_old
  }
  worst_value <- max(good_scalers)
  diff_hist <- list(good_points = good_points, good_values = good_values,
      failed_points = failed_points)
  svar <- worst_value - best_value
  sflag <- 1
  if (abs(best_value_old - best_value) < 1e-14) {
    sflag <- 0
    jac <- c()
    grad <- c()
  }
  list(sflag = sflag, best_value = best_value, best_value_f = best_value_f,
  best_point = best_point, svar = svar, diff_hist = diff_hist)
}