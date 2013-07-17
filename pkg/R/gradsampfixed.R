gradsampfixed <-
function(fn, gr, x0, nvar=length(x0), samprad, f0=fn(x0), g0=gr(x0), maxit=1000,normtol=1e-6, ngrad=min(2*nvar,100,nvar+10), fvalquit=-Inf, prtlevel=1){
  
  #initialisations and declarations
  x <- x0
  f <- f0
  g <- g0
  X <- x
  G <- g
  w <- 1
  quitall <- FALSE
  dnorm <- Inf
  
  for(iter in 1:maxit){
    tmp <- getbundle(fn, gr, x, g, samprad, ngrad)
    Xnew <- tmp$xbundle
    Gnew <- tmp$gbundle
    tmp <- qpspecial(Gnew)
    wnew <- tmp$x
    dnew <- tmp$d
    
    dnew <- -dnew
    gtdnew <- sum(g*dnew)
    dnormnew <- sqrt(sum(dnew*dnew))
    if(dnormnew<dnorm){
      dnorm <- dnormnew
      X <- Xnew
      G <- Gnew
      w <- wnew
    }
    if(dnormnew<normtol){
      mess <- paste("gradsamp: tolerance met at iter=",iter,"\n")
      if(prtlevel) cat(mess)
      return(list(x=x, f=f, g=g, dnorm=dnorm, X=X, G=G, w=w, quitall=quitall,
		  message = mess))
    }
    else if(gtdnew >=0 || is.na(gtdnew)){
	   mess <- paste("gradsamp: not descent direction, quit at iter=",iter,"\n")
     if(prtlevel) cat(mess)
       return(list(x=x, f=f, g=g, dnorm=dnorm, X=X, G=G, w=w, quitall=quitall,
		  message = mess))
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
      quitall <-TRUE
      mess <- paste("gradsamp: reached target objective, quit at iter=",iter,"\n")
      if(prtlevel) cat(mess)
       return(list(x=x, f=f, g=g, dnorm=dnorm, X=X, G=G, w=w, quitall=quitall,
		  message = mess))
    }
    if(fail == -1){
      mess <- paste("gradsamp: f may be unbounded below, quit at iter =",iter,"\n")
      if(prtlevel) cat(mess)
      quitall=TRUE
       return(list(x=x, f=f, g=g, dnorm=dnorm, X=X, G=G, w=w, quitall=quitall,
		  message = mess))
    }
    #ommit cpu conditions
  }
  mess <- paste("gradsamp: completed iteration=",iter,"\n")  
  if(prtlevel) cat(mess)
   return(list(x=x, f=f, g=g, dnorm=dnorm, X=X, G=G, w=w, quitall=quitall,
		  message = mess))
}
