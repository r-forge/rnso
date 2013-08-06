imfil <- function(x0, fn, budget, bounds, options = imfil_optset()){
  qbounds <- bounds
  dbounds <- bounds[, 2] - bounds[, 1]
  n <- length(x0)
  imfil_smooth_problem <- options$smooth_problem
  imfil_complete_history <- options$_complete_history
  if (imfil_smooth_problem){
    #TODO
  }
  imfil_fscale <- options$fscale
  
}