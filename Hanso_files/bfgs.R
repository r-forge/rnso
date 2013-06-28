bfgs <- function(fn,gr,nvar,maxit = 1000, normtol = 1e-6, 
		      fvalquit = -Inf, xnormquit = Inf, nvec = 0, prtlevel = 1,
		      strongwolfe = 0, wolfe1 = 1e-4, wolfe2 = 0.5, quitLSfail = 1,
		      ngrad = 2, evaldist = 1e-4, H0 = diag(nvar), scale = 1,nstart=10,
		      x0 = matrix(rnorm(nvar*nstart),nvar,nstart){
  x <- list()
  f <- c()
  d <- list()
  iter <- c()
  #HH <- list()
  H <- list()
  info <- c()
  X <- list()
  G <- list()
  w <- list()
  fevalrec <- c()
  xrec <- list()
  Hrec <- list()
  
  for(run in 1:nstart){
    tmp <- bfgs1run(fn,gr,options$x0[,run],options)
    #print(tmp)
    x[[run]] <- tmp$x
    f[run] <- tmp$f
    d[[run]] <- tmp$d
    HH <- as.matrix(tmp$H)
    iter[run] <- tmp$iter
    info[run] <- tmp$info
    X[[run]] <- tmp$X
    G[[run]] <- tmp$G
    w[[run]] <- tmp$w
    fevalrec[run] <- tmp$fevalrec
    xrec[[run]] <- tmp$xrec
    Hrec[[run]] <- tmp$Hrec
  }
  H[[run]] <- (HH+t(HH))/2
  #cputime break condition
  if(nstart == 1){
    H <- H[1]
    fevalrec <- fevalrec[1]
    xrec <- xrec[1]
    Hrec <- Hrec[1]
    X <- X[1]
    G <- G[1]
    w <- w[1]
  }
  return(list(x,f,d,H,iter,X,G,fevalrec,xrec,Hrec))
}