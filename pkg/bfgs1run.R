bfgs1run=function(x0,pars,options){
  n=pars$nvar
  fgname=pars$fgname
  normtol=options$normtol
  fvalquit=options$fvalquit
  xnormquit=options$xnormquit
  cpufinish=cputime+options$cpumax
  maxit=options$maxit
  nvec=options$nvec
  prtlevel=options$prtlevel
  strongwolfe=options$strongwolfe
  wolfe1=options$wolfe1
  wolfe2=options$wolfe2
  quitLSfail=options$quitLSfail
  ngrad=options$ngrad
  evaldist=options$evaldist
  H0=options$H0
  H=H0
  scale=options$scale
  x=x0
  tmp=eval(fgname(x,pars))
  f=tmp$f
  g=tmp$g
  d=g
  G=g
  X=x
  nG=1
  w=1
  dnorm=norm(g,type="1")
  if(nvec >0){
    S=c()
    Y=c()    
  }
  iter=0
  ##print error msgs
  for(iter=1:maxit){
    if(nvec==0){
      p=-H*g
    }else
      p=-hgprod(H,g,S,Y)
    gtp=t(g)%*%p
    if(gtp>=0 | is.nan(gtp)){
      if(prtlevel >0) print("bfgs non descent...")
    }info=6
    return()
  }
  gprev=g
  if(strongwolfe){
    fevalrecline=nan;
    tmp=linesch_sw(x,f,g,p,pars,wolfe1,wolfe2,fvalquit,prtlevel)
    if(wolfe2==0){
      increase=1e-8*(1+alpha)
      x=x+increase*p
    }
    if(prtlevel>1) print("exact line sch simulation slightly increasing....")
    tmp=eval(pars.fgname,x,pars)
  }
  else {
    tmp=linesch_ww(x,f,g,p,pars,wolfe1,wolfe2,fvalquit,prtlevel)
  }
  if(alpha*norm(p)>evaldist){
    nG=1
    G=g
    X=x
  }
  if(nG<ngrad){
    G=rbind(g,G[,1:ngrad-1])
    X=rbind(x,X[,1:ngrad-1])
  }
  if(nG>1) {
    qtemp=qpspecial(G)
    w=qtemp$w
    d=qtemp$d
  }else{
   w=1
   d=g
  }
  dnorm=norm(d)
  xrec[,iter]=x
  fevalrec[iter]=fevalrecline
  Hrec[iter]=H
  if(prtlevel>1){
    nfeval=length(fevalrecline)
    print("bfgs:...")
  }
  if(f<fvalquit){
    if(prtlevel>0) print("bfgs reached objective")
  }
  info=2
  return(?)
  if(norm(X)>xnormquit){
    if(prtlevel>0) print("bfgs: norm exceeds specified limit")
    info=3
    return(?)
  }
  if(fail==1){
    if(~quitLSfail){
      if(prtlevel>1) print("bfgs: continue alhough line search failed")
    }else{
      if(prtlevel>0) print("bfgs: quit at iteration...")
      info=7
      return(?)
    }
  }
  if(fail==-1){
    if(prtlevel>0) print("bfgs: f may be unbounded below...")
    info=8
    return(?)
  }
  if(dnorm<=normtol){
    if(prtlevel>0){
      if(nG==1) print("bfgs: gradient norm below tolerance, quite iteration")
      else print("bfgs: norm of smallest vector in convex hull of gradients below tolerance, quit...")
    }
    info=0
    return(?)
  }
  if(cputime>cpufinish){
    if(prtlevel>0) print("bfgs: cpu time limit exceeded, quit iteration")
  
    info=4
    return(?)
  }
  s=alpha*p
  y=g-gprev
  sty=t(s)%*%y
  if(nvec==0){
    if(sty>0){
      if(iter==1 & scale) H=(sty/(t(y)%*%y))*H
      rho=1/sty
      rhoHyst=rho*(H*y)*t(s)
      H=H-t(rhoHyst)-rhoHyst+rho*s*(t(y)%*%rhoHyst)+rho*s%*%t(s)
    }else{
      if(prtlevel>1) print("bfgs: sty<=0, skipping bfgs update at iteration...")
    }
  }else{
    s=alpha*p
    y=g-gprev
    if(iter<=nvec){
      S=rbind(S[,2:nvec],s)
      Y=rbind(Y[,2:nvec],y)
    }
    if(scale) H=((t(s)%*%y)/(t(y)%*%y))*H0
  }
  
}





















