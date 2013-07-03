bfgs <- function(fn,gr,nvar,nstart=10,x0 = matrix(rnorm(nvar*nstart),nvar,nstart),maxit = 1000, normtol = 1e-6, 
		      fvalquit = -Inf, xnormquit = Inf, nvec = 0, prtlevel = 1,
		      strongwolfe = 0, wolfe1 = 1e-4, wolfe2 = 0.5, quitLSfail = 1,
		      ngrad = 0, evaldist = 1e-4, H0 = NULL, scale = 1)
		      {
  x <- matrix(NA,nvar,nstart)
  f <- c()
  d <- list()
  iter <- c()
  
  HH <-list()
  info <- c()
  X <- list()
  G <- list()
  w <- list()
  mess <- c()
  
  for(run in 1:nstart){
    tmp <- bfgs1run(fn,gr,x0[,run],H0, maxit,  fvalquit,
                    normtol, xnormquit, evaldist, ngrad,
                    scale, wolfe1, wolfe2, quitLSfail)
    #print(tmp)
    x[,run] <- tmp$x
    f[run] <- tmp$f
    d[[run]] <- tmp$d
    H <- as.matrix(tmp$H)
    iter[run] <- tmp$iter
    info[run] <- tmp$info
    mess[run] <- tmp$message
    X[[run]] <- tmp$X
    G[[run]] <- tmp$G
    w[[run]] <- tmp$w
    HH[[run]] <- (H+t(H))/2
  }
  #no need for cpu break
  #no need for special formatting for nstart==1
#   if(nstart == 1){
#     H <- H[1]
#     fevalrec <- fevalrec[1]
#     xrec <- xrec[1]
#     Hrec <- Hrec[1]
#     X <- X[1]
#     G <- G[1]
#     w <- w[1]
#   }
  return(list(x=x,f=f,d=d,H=HH,iter=iter,message=mess,X=X,G=G,w=w))
}