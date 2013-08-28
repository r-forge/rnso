imfil_create_stencil <- function(options, n){
  stencil  <- options$stencil
  vstencil <- options$vstencil
  if (!is.null(vstencil)){
    stencil <- -1
  }
  if (stencil == -1){
    
  } else if (stencil == 0) {
    v <- cbind(diag(n),-diag(n))
  } else if (stencil == 1) {
    v <- diag(n)
  } else if (stencil == 2) {
    v <- cbind(diag(n), -matrix(1,n,1)/sqrt(n))
  } else stop("imfil_create_stencil: illegal stencil")
  if (stencil != -1){
    vstencil <- v
  }
  vstencil
}
