gradsampfixed <- function(fn, gr, x0, nvar=length(x0), samprad, f0=fn(x0), g0=gr(x0), maxit=1000,normtol=1e-6, ngrad=min(2*nvar,100,nvar+10), fvalquit=-Inf, prtlevel=1){
  
  #initialisations and declarations
  x <- x0
  f <- f0
  g <- g0
  X <- x
  G <- g
  w <- 1
  quitall <- 0
  dnorm <- Inf
  
  for(iter in 1:maxit){
    tmp <- getbundle(fn, gr, x, g, samprad, ngrad)
    Xnew <- tmp$xbundle
    Gnew <- tmp$gbundle
    
    tmp <- qpspecial(Gnew)
    wnew <- tmp$x
    dnew <- tmp$d
    
    dnew <- -dnew
    gtdnew <- t(g)%*%dnew
    dnormnew <- norm(dnew)
    if(dnormnew<dnorm){
      dnorm <- dnormnew
      X <- Xnew
      G <- Gnew
      w <- wnew
    }
    if(dnormnew<normtol){
      return(list(x=x, f=f, g=g, dnorm=dnorm, X=X, G=G, w=w, quitall=quitall))
    }
    else if(gtdnew >=0 | is.nan(gtdnew)){
      if (prtlevel>0) warning("not descent direction, quit")
      return(list(x=x, f=f, g=g, dnorm=dnorm, X=X, G=G, w=w, quitall=quitall))
    }
    wolfe1 <- 0
    wolfe2 <- 0
   
    tmp <- linesch_ww(fn,gr,x,dnew,f,g,wolfe1,wolfe2,fvalquit,prtlevel)
    alpha <- tmp$alpha
    x <- tmp$xalpha
    f <- tmp$falpha
    g <- tmp$galpha
    fail <- tmp$fail
    
    if(f< fvalquit){
      quitall <-1
      return(list(x=x, f=f, g=g, dnorm=dnorm, X=X, G=G, w=w, quitall=quitall))
    }
    if(fail == -1){
      quitall=1
      return(list(x=x, f=f, g=g, dnorm=dnorm, X=X, G=G, w=w, quitall=quitall))
    }
    #ommit cpu conditions
  }
  return(list(x=x, f=f, g=g, dnorm=dnorm, X=X, G=G, w=w, quitall=quitall))
}  
    
