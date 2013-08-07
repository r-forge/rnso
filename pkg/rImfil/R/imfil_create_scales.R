imfil_create_scales <- function(options){
  custom_scales <- options$custom_scales
  mcs <- length(custom_scales)
  if (mcs > 0){
    dscal <- custom_scales
  } else {
    scalestart <- options$scalestart
    scaledepth <- options$scaledepth
  }
  if (scalestart > scaledepth) stop("imfil_create_scales: scalestart > scaledepth")
  dscal <- -(scalestart : scaledepth)
  dscal <- 2^dscal
  dscal
}