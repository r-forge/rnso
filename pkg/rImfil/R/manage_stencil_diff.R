manage_stencil_diff <- function(x,fn,funs,iteration_data,fcount, stencil_data,stop_now) {
  
  if (stop_now == 0) {
    options <- iteration_data$options
    imfil_complete_history <- options$complete_history
    h <- iteration_data$h
    stencil_delta <- stencil_data$stencil_delta
    svarmin <- stencil_data$svarmin
    noise_val <- stencil_data$noise_val
    obounds <- iteration_data$bounds
    bounds <- stencil_data$bounds
    v <- stencil_data$v
    vv <- imfil_augment_directions(x, v, h, options, bounds)
    complete_history <- iteration_data$complete_history
    tmp <- stencil_diff(x, fn, h*vv, funs, iteration_data, complete_history)
    sgrad <- tmp$grad
    fb <- tmp$best_value
    fbf <- tmp$best_value_f
    xb <- tmp$best_point
    icount <- tmp$icount
    sflag <- tmp$sflag
    svar <- tmp$svar
    diff_hist <- tmp$diff_hist
    jac <- tmp$jac
    
    focunt <- fcount + icount
    pgrad <- x - kk_proj(x - sgrad, obounds)
    npgrad <- norm(pgrad, "I")
    #browser()
    if(is.na(svar)) svar<- 0 # fix for svar being NA
    if (!is.na(svar)){ if (max(noise_val, svarmin) > svar) sflag <- 0} else{
      sflag <- 0
    } #fix
    if (!is.na(svar)){ if (stencil_delta > svar) {
      stop_now <- 1
      sflag <- 0
    }} else{
      stop_now <- 1
      sflag <- 0
    }
    iteration_datap <- iteration_data
    if (imfil_complete_history > 0) {
      complete_history <- many_point_hist_update(complete_history, diff_hist)
      iteration_datap$complete_history <- complete_history
    }
    tmp <- reconcile_best_point(fbf, xb, iteration_datap)
    iteration_datap <- tmp$new_data
    rflag <- tmp$rflag
    
    sdiff <- jac_or_grad(sgrad, jac, options)
    
    list(sdiff = sdiff, sgrad = sgrad, npgrad = npgrad, fcount = fcount, 
    sflag = sflag, jac = jac, iteration_datap = iteration_datap, stop_now = stop_now)
    }
}