imfil <- function(x0, fn, budget, bounds, options = imfil_optset()){
  qbounds <- bounds
  dbounds <- bounds[, 2] - bounds[, 1]
  n <- length(x0)
  imfil_smooth_problem <- options$smooth_problem
  imfil_complete_history <- options$complete_history
  if (imfil_smooth_problem == 1){
    #TODO
  }
  imfil_fscale <- options$fscale
  fun_data <- list(fn = fn, qbounds = qbounds, dbounds = dbounds, 
		    fun_scale = imfil_fscale)
  core_data <- list(options = options, fun_data = fun_data)
  
  #badargs <- imfil_error_check() #TODO
  
  #if(badargs)
  
  z0 <- (x0 -qbounds[,1])/dbounds
  tmp <- imfil_core(z0, f_internal, budget, core_data, bounds)
  z <- tmp$
  histout <- tmp$
  complete_history <- tmp$
  ifail <- tmp$
}