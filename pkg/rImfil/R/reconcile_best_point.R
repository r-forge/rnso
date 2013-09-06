reconcile_best_point <- function(funs, x, old_data) {
  new_data <- old_data
  least_squares <- old_data$options$least_squares
  rflag <- 1
  fb <- old_data$fobjb
  fval <- f_to_vals(funs, least_squares)
  if (fval < fb) {
    rflag <- 0
    new_data$xb <- x
    new_data$funsb <- funs
    new_data$fobjb <- fval
  }
  list(new_data = new_data, rflag = rflag)
}