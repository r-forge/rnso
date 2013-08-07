stencil_diff <- function(x, fn, dx, fc, iteration_data, complete_history) {
  core_data <- interation_data$core_data
  h <- iteration_data$h
  bounds <- iteration_data$bounds
  options <- core_data$options
  parallel <- options$parallel
  least_squares <- options$least_squares
  tmp <- imfil_poll_stencil(x, fn, dx, fc, bounds, core_data, h, complete_history)
  best_value <- tmp$oldindex
}