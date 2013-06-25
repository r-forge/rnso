gradsamp1run=function(fn,gr,x0,f0=fn(x0), g0=gr(x0), samprad=c(1e-3,1e-4,1e-4),maxit=100,normtol=1e-4, ngrad=2, fvalquit=-Inf, prtlevel=1){

  for(choice in 1:length(samprad)){
    tmp=gradsampfixed(fn,gr,x0,f0,g0,samprad[choice],maxit,normtol,ngrad,fvalquit,prtlevel)
    #print(tmp)
    if (tmp$quitall) return(list(x=tmp$x,f=tmp$f,g=tmp$g,dnorm=tmp$dnorm,X=tmp$X,G=tmp$G,w=tmp$w))
    x0=tmp$x
   
    f0=tmp$f
    g0=tmp$g
  }
  return(list(x=x0,f=f0,g=g0,dnorm=tmp$dnorm,X=tmp$X,G=tmp$G,w=tmp$w))
}