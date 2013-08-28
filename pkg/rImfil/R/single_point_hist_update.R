single_point_hist_update <- function(old_hist, x, fout, ifail){
  new_hist <- old_hist
  if (ifail == 1){
    if (is.na(new_hist$failed_points)) {
      new_hist$failed_points <- old_hist$failed_points
    } else {
      new_hist$failed_points <- cbind(old_hist$failed_points, x)
    }
  } else {
    if (is.na(new_hist$good_points)) {
        new_hist$good_points <- old_hist$good_points
      }else {
        new_hist$good_points   <- cbind(old_hist$good_points, x)
      }
    if (is.na(new_hist$good_values)){
      new_hist$good_values <- old_hist$good_values
    }else{
      new_hist$good_values   <- cbind(old_hist$good_values, fout)
    }
  }
  new_hist
}