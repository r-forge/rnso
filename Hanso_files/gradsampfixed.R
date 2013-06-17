gradsampfixed=function(x0,f0,g0,samprad,pars,options){
  prtlevel=options$prtlevel
  x=x0
  f=f0
  g=g0
  X=x
  G=g
  w=1
  quitall=0
  maxit=options$maxit
  normtol=options$normtol
  ngrad=options$ngrad
  fvalquit=options$fvalquit
  #
  dnorm=Inf
  for(iter=1:maxit){
    tmp=getbundle(x,g,samprad,ngrad,pars)
    Xnew=tmp$xbundle
    Gnew=tmp$gbundle
    tmp=qpspecial(Gnew)
    wnew=tmp$x
    dnew=tmp$d
    dnew=-dnew
    gtdnew=t(g)%*%dnew
    dnormnew=norm(dnew)
    if(dnormnew<dnorm){
      dnorm=dnormnew
      X=Xnew
      G=Gnew
      w=wnew
    }
    if(dnormnew<normtol){
      #return? with warnings
    }
    wolfe1=0
    wolfe2=0
    tmp=linesch_ww(x,f,g,dnew,pars,wolfe1,wolfe2,fvalquit,prtlevel)
    alpha=tmp$alpha
    x=tmp$xalpha
    f=tmp$falpha
    g=tmp$gradalpha
    fail=tmp$fail
    #quitall conditions
    #warnings
    
    
    
    
    
    
    
    
    
    
  }
  return(list(x,g,dnorm,X,G,w,quitall))
}
