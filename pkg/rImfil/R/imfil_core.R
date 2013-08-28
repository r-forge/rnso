imfil_core <- function (x0, fn , budget, core_data, bounds) {
  #global imfil_fscale
  fcount <- 0
  ifailed <- 0
  options <- core_data$options
  n <- length(x0)
  obounds <- matrix(0, n, 2)
  obounds[, 2] <- 1
  
  imfil_explore <- options$explore
  imfil_executive <- options$executive
  imfil_function_delta <- options$function_delta
  imfil_least_squares <- options$least_squares
  imfil_maxit <- options$maxit
  imfil_maxitarm <- options$maxitarm
  imfil_maxfail <- options$maxfail
  imfil_noise_aware <- options$noise_aware
  imfil_scale_aware <- options$scale_aware
  imfil_target <- options$target
  imfil_termtol <- options$termtol
  stencil_wins <- options$stencil_wins
  verbose <- options$verbose
  imfil_complete_history <- options$complete_history
  complete_history <- list(good_points = matrix(), good_values = matrix(), 
      failed_points = matrix())
  
  tmp <- setup_explore(options)
  explore_function <- tmp$explore_function
  explore_data_flag <- tmp$explore_data_flag
  explore_data <- tmp$explore_data
  x <- x0
  xold <- x0
  n <- length(x0)
  histout <- matrix()
  dscal <- imfil_create_scales(options)
  nscal <- length(dscal)
  imfil_exec <- options$executive
  
  if (imfil_exec == 1) {
    hess <- options$executive_data
  } else {
    hess <- diag(n)
  }
  xc <- x0
  ns <- 0
  failc <- 0
  stop_now <- 0
  fval <- imfil_target + 1
  
  sflag <- 1
  while (ns < nscal && fcount <= budget && imfil_maxfail && stop_now == 0 && fval > imfil_target) {
    ns <- ns + 1
    h <- dscal[ns]
    if (imfil_noise_aware > 0 || fcount == 0 || imfil_scale_aware > 0) {
      tmp <- f_internal(x, h, core_data)
      funs <- tmp$fx
      iff <- tmp$iff
      icf <- tmp$icf
      noise_val <- tmp$tol # remember! fn is f_internal not f_easy
      fval <- f_to_vals(funs, imfil_least_squares)
      icount <- icf
    } else {
      icount <- 0
    }
    
    if (fcount == 0) {
      #funerr <- imfil_error_check() 
      imfil_target <- imfil_target/imfil_fscale
      imfil_function_delta <- imfil_function_delta/imfil_fscale
      if (is.na(histout)){
        histout <- rbind(1, fval, 0, 0, 0, as.matrix(x))
      } else {
        histout <- rbind(as.numeric(histout),1, fval, 0, 0, 0, as.matrix(x))
      }
      
      
      itc <- 0
      stencil_data <- create_stencil_data(options, imfil_fscale, noise_val, bounds)
      iteration_data <- list(h = h, obounds = obounds, itc = itc, xb = x,
	fobjb = fval, funsb = funs, complete_history = complete_history, 
	f_internal = fn, core_data = core_data, options = options)
    } else {
      iteration_data$h <- h
      stencil_data$noise_val <- noise_val
    }
    if (imfil_complete_history > 0) {
      complete_history <- iteration_data$complete_history
      complete_history <- single_point_hist_update(complete_history, x, funs, iff)
      iteration_data$complete_history <- complete_history
    }
    fcount <- fcount + icf
    if (fval < imfil_target) {
      stop_now <- 1
      break
    }
    stol <- imfil_termtol*h
    iarm <- 0
    nfail <- 0
    
    tmp <- manage_stencil_diff(x, fn, funs, iteration_data, fcount, stencil_data, stop_now)
    sdiff <- tmp$sdiff
    sgrad <- tmp$sgrad
    npgrad <- tmp$npgrad
    fcount <- tmp$fcount
    sflag <- tmp$sflag
    jac <- tmp$jac
    iteration_data <- tmp$iteration_datap
    stop_now <- tmp$stop_now
    histout <- t(cbind(t(histout), t(cbind(fcount, fval, npgrad, -1, -1, t(x)))))
    gc <- sgrad
  }
  
  if (npgrad < stol || sflag == 0 || stop_now == 1) {
    gc <- sgrad
    if (sflag != 0) {
      failc <- failc + 1
    } else {
      failc <- 0
    }
  } else {
    failc <- 0
    itc <- 0
    while (itc < imfil_maxit*n && fval > imfil_target && npgrad >= stol && nfail == 0 && fcount < budget && sflag > 0) {
      itc <- itc + 1
      iteration_data$itc <- itc
      fc <- fval
      tmp <- imfil_inner_iteration (fn, x, funs, sdiff, xc, gc, iteration_data, hess, fcount)
      xp <- tmp$xp
      fval <- tmp$fval
      funs <- tmp$funs
      fcount <- tmp$fcount
      iarm <- tmp$iarm
      iteration_data <- tmp$iteration_datap
      nfail <- tmp$nfail
      hess <- tmp$hess
      
      if (fval < imfil_target) {
	stop_now <- 1
	x <- xp
	stepn <- norm(xold-x, "I")
	histout <- t(cbind(t(histout), t(cbind(fcount, fval, npgrad, stepn, iarm, t(x)))))
	break
      }
      xold <- x
      xc <- x
      gc <- sgrad
      x <- xp
      if (stencil_wins == 1) {
	
      }
    }
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  list(x = x, histout = histout, complete_history = complete_history, ifailed = ifailed)
}