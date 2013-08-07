single_point_hist_update <- function(old_hist, x, fout, ifail){
  new_hist <- old_hist
  if (ifail == 1){
    new_hist$failed_points <- cbind(old_hist$failed_points, x)
  } else {
    new_hist$good_points   <- cbind(old_hist$good_points, x)
    new_hist$good_values   <- cbind(old_hist$good_values, fout)
  }
  new_hist
}