write_best_to_x <- function(iteration_data) {
 x <- iteration_data$xb
 funs <- iteration_data$funsb
 fval <- iteration_data$fobjb
 list(x = x, funs = funs, fval = fval)
 }