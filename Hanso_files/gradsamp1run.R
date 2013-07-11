gradsamp1run=function(fn,gr,x0,nvar=length(x0),f0=fn(x0), g0=gr(x0), samprad=c(1e-4,1e-5,1e-6),maxit=1000,normtol=1e-4, ngrad=min(2*nvar,100,nvar+10), fvalquit=-Inf, prtlevel=1){

  mess <- c()
  for(choice in 1:length(samprad)){
    tmp <- gradsampfixed(fn,gr,x0,nvar,samprad[choice],f0,g0,maxit,normtol,ngrad,fvalquit,prtlevel)
    #print(tmp)
    if (tmp$quitall) return(list(x=tmp$x,f=tmp$f,g=tmp$g,dnorm=tmp$dnorm,X=tmp$X,G=tmp$G,w=tmp$w))
    x0 <- tmp$x
   
    f0 <- tmp$f
    g0 <- tmp$g
    mess[choice] <- tmp$mess
  }
  return(list(x=x0,f=f0,g=g0,dnorm=tmp$dnorm,X=tmp$X,G=tmp$G,w=tmp$w,
	    message = c(mess)))
}