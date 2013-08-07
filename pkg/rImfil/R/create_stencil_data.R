create_stencil_data <- function(options, imfil_fscale, noise_val, bounds){
  n <- length(bounds[, 1])
  v <- imfil_create_stencil(options, n)
  imfil_stencil_delta <- options$stencil_delta
  imfil_svarmin       <- options$svarmin
  imfil_stencil_delta <- imfil_stencil_delta/imfil_fscale
  imfil_svarmin       <- imfil_svarmin/imfil_fscale
  
  stencil_data <- list(stencil_delta = imfil_stencil_delta, 
			svarmin = imfil_svarmin, 
			noise_val = noise_val, 
			bounds = bounds)
  stencil_data			
}