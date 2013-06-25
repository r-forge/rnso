getbundle <- function(fn,gr,x,g=gr(x),samprad,N){
  m <- length(x)
  #declare empty matrices
  xbundle <- matrix(NA,m,N)
  gbundle <- matrix(NA,m,N)
  xbundle[,1] <- x
  gbundle[,1] <- g
  for(k in 2:N){
    xpert <- x+samprad*(runif(m)-0.5) #samprad is a scaler here
    f <- fn(x)
    grad <- gr(x)
    count <- 0
    while(isnaninf(f) | isnaninf(grad)){
      xpert <- (x+xpert)/2
      f <- fn(x)
      grad <- gr(x)
     # count <- count+1
    }
    xbundle[,k] <- xpert
    gbundle[,k] <- grad
  }
  return(list(xbundle=xbundle,gbundle=gbundle))
}