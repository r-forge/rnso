many_point_hist_update <- function(old_hist, diff_hist) {
  new_hist <- old_hist
  #if (!is.na(new_hist$good_points)){
  #  new_hist$good_points   <- cbind(new_hist$good_points, diff_hist$good_points)
  #} else {
  new_hist$good_points <- diff_hist$good_points #}
  #if (!is.na(new_hist$good_values)){
   # new_hist$good_values   <- cbind(new_hist$good_values, diff_hist$good_values)
  #} else {
  new_hist$good_values <- diff_hist$good_values #}
  #if (!is.na(new_hist$failed_points)){
   # new_hist$failed_points   <- cbind(new_hist$failed_points, diff_hist$failed_points)
  #} else {
  new_hist$failed_points <- diff_hist$failed_points #}
  
  new_hist
}