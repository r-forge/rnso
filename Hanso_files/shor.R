shor <- function(fn,gr,x0, maxit = 1000,  fvalquit = -Inf,
                     normtol = 1e-6, xnormquit = Inf, evaldist = 1e-4, 
                     ngrad = 0, rescale = 0, strongwolfe = 0, useprevstep = 0,
                     wolfe1 = 1e-4, wolfe2 = 0.5, quitLSfail = TRUE,prtlevel=1)
                     {

  x <- matrix(NA,nvar,nstart)
  f <- c()
  g <- matrix(NA,nvar,nstart)
  B <- list()
  frec <- list()
  fevalrec <- list()
  betarec <- list()
  xrec <- list()
  svrec <- list()
  
  for(run in 1:nstart){
    res <- shor1run()
    x[,run] <- res$x
    f[run] <- res$f
    g[,run] <-res$g
    B[[run]] <- res$B
    frec[[run]] <- res$frec
    fevalrec[[run]] <- res$fevalrecall
    betarec[[run]] <- res$betarec
    xrec[[run]] <- res$xrec
    svrec[[run]] <- res$svrec
    
  }
  
  return(list(x = x, f = f, g = g, B = B, frec = frec, fevalrec = fevalrec,
  betarec = betarec, xrec = xrec, svrec = svrec))
  
}