imfil_augment_directions <- function(x, vin, h, options, bounds) {
  random_stencil <- options$random_stencil
  vout <- random_augment(vin, random_stencil)
  new_directions <- options$add_new_directions
  lnew <- length(new_directions)
  
  if (lnew > 0) {
    dbv <- bounds[, 2] - bounds[, 1]
    db <- diag(dbv)
    unscaled_x <- db%*%x + bounds[, 1]
    unscaled_v <- db%*%vout
    unscaled_vnew <- do.call(new_directions, unscaled_x, h, unscaled_v)
    mv <- nrow(unscaled_vnew)
    nv <- ncol(unscaled_vnew)
    if (mv > 0) {
      vnew <- solve(db)%*%unscaled_vnew
      mv <- nrow(vnew)
      nv <- ncol(vnew)
      for (i in 1:nv){
	vnew[, i] <- vnew[, i]/norm(vnew[, i],"I")
      }
      vout <- cbind(vout, vnew)
    }
  }
  vout
}