manage_explore <- function (explore_function,fn,iteration_data,fctin,explore_data){
  iteration_datap <- iteration_data
  options <- iteration_data$options
  imfil_complete_history <- options(complete_history)
  tmp <- explore_function(fn, iteration_data, explore_data)
  xs <- tmp$x
  fs <- tmp$fv
  my_cost <- tmp$cost
  explore_hist <- tmp$hist #check explore function interface
  fctout <- fctin + my_cost
  
  if (imfil_complete_history > 0) {
    iteration_datap$complete_history <- many_point_hist_update(iteration_datap$complete_history,
     explore_hist)
  }
  tmp <- reconcile_best_point(fs, xs, iteration_datap)
  iteration_datap <- tmp$new_data
  rflag <- tmp$rflag
  list(fctout = fctout, iteration_datap = iteration_datap)
}