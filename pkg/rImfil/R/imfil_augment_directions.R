imfil_augment_directions <- function(x, vin, h, options, bounds) {
  random_stencil <- options$random_stencil
  vout <- random_augment(vin, random_stencil)
  new_directions <- options$add_new_directions
  lnew <- length(new_directions)
}