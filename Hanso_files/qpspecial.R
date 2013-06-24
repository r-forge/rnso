qpspecial <- function(G,x = matrix(1,ncol(G),1),maxit = 100){
  m <- nrow(G)
  n <- ncol(G)
  if(length(G)  == 0){
    
    Warning('G is empty')
    Return(list(x <- c(), q <- Inf, Info <- c(2,0))) 
  }
  
  maxit <- max(maxit,10) 

  e <- matrix(1,n,1)
  x <- as.matrix(c(x))
  nx <- length(x)
  if(any(x < 0) || nx !=  n) {x <- e}
  
  idx <- seq(1,(n*n),by <- n+1)
  Q <- t(G)%*%G
  z <- x
  y <- 0
  eta <- 0.9995
  delta <- 3
  mu0 <- sum(x*z)/n
  
  tolmu <- 1e-5
  tolrs <- 1e-5
  kmu <- tolmu*mu0
  
  nQ <- norm(Q,"I")+2
  krs <- tolrs*nQ
  ap <- 0
  ad <- 0
  for(k in 1:maxit){
    r1 <- -Q%*%x+e%*%y+z
    r2 <- -1+sum(x)
    r3 <- -x*z
    rs <- norm(rbind(r1,r2),"I") #check norm(,inf)?
    
    mu <- -sum(r3)/n
    #cat("iteration= ",k ,"\t mu= ",mu,"\t rs= ",rs,"\t krs= ",krs,"\n") #for debugging
    if(mu < kmu){
      if(rs<krs) {info <- c(0,k-1); break}
    }
    zdx <- z/x
    QD <- Q
    QD[idx] <- QD[idx]+zdx
    C <- chol(QD)
    KT <- solve(t(C) , e)
    M <- t(KT)%*%KT
    r4 <- r1+r3/x
    r5 <- t(KT)%*%(solve(t(C),r4))
    r6 <- r2+r5
    dy <- -r6/M
    r7 <- r4+e%*%dy
    dx <- solve(C,(solve(t(C),r7)))
    dz <- (r3-z*dx)/x
    p <- -x/dx
    ap <- min(min(p[p>0]),1)
    if(length(ap) == 0) ap <- 1
    p <- -z/dz
    ad <- min(min(p[p>0]),1)
    if(length(ad) == 0) ad <- 1
    mauff <- (t(x+ap*dx)%*%(z+ad*dz))/n
    sig <- (mauff/mu)^delta
    r3 <- r3+rep(sig*mu,n)
    r3 <- r3-dx*dz
    r4 <- r1+r3/x
    r5 <- t(KT)%*%(solve(t(C),r4))
    r6 <- r2+r5
    dy <- -r6/M
    r7 <- r4+e%*%dy
    dx <- solve(C,(solve(t(C),r7)))
    dz <- (r3-z*dx)/x
    p <- -x/dx
    ap <- min(min(p[p>0]),1)
    if(length(ap) == 0) ap <- 1
    p <- -z/dz
    ad <- min(min(p[p>0]),1)
    if(length(ad) == 0) ad <- 1
    x <- x+eta*ap*dx
    y <- y+eta*ad*dy
    z <- z+eta*ad*dz
  }
  if(k >= maxit) info <- 1
  x <- pmax(x,0)
  x <- x/sum(x)
  d <- G%*%x
  q <- t(d)%*%d
  return(list(x <- x,d <- d,q <- q,info <- info))
  
}























