bfgs1run=function(fn,gr,x0,options){
  n=length(x0)
  normtol=options$normtol
  fvalquit=options$fvalquit
  xnormquit=options$xnormquit
  #cpufinish=cputime+options$cpumax
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
  f=fn(x)
  g=gr(x)
  d=g
  G=g
  X=x
  nG=1
  w=1
  dnorm=norm(as.matrix(g),type="1")
  if(nvec >0){
    S=c()
    Y=c()    
  }
  iter=0
  #initializations
  fevalrec=list()
  xrec=list()
  Hrec=list()
  ##print error msgs
  for(iter in 1:maxit){
    if(nvec==0){
      p=-H%*%g
    }else
      p=-hgprod(H,g,S,Y) #why!! it i*s just computing H%*%g 
    gtp=t(g)%*%p
    if(gtp>=0 | is.nan(gtp)){
      if(prtlevel >0) stop("bfgs non descent direction, quit")
      info=6
      #      return(list(x=x,f=f,d=d,H=H,iter=iter,info=info,X=X,G=G,w=w,fevalrec=fevalrec,xrec=xrec,Hrec=Hrec))
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
      tmp=linesch_ww(fn,gr,x0,p)
      alpha=tmp$alpha
      x=tmp$xalpha
      f=tmp$falpha
      g=tmp$galpha
      fail=tmp$fail
      fevalrecline=tmp$fevalrec
    }
    if(alpha*norm(p)>evaldist){
      nG=1
      G=g
      X=x
    }
    else if(nG<ngrad){
      #something wrong here
      nG=nG+1
      G=cbind(g,G)
      X=cbind(x,X)
    }
    else{
      G=cbind(g,G[,1:ngrad-1])
      X=cbind(x,X[,1:ngrad-1])
    }
    
    
    if(nG>1) {
      qtemp=qpspecial(t(G))
      w=qtemp$x
      d=qtemp$d
    }else{
      w=1
      d=g
    }
    dnorm=norm(as.matrix(d))
    xrec[[iter]]=x #xrec is a list
    fevalrec[[iter]]=fevalrecline
    Hrec[[iter]]=H
    if(prtlevel>1){
      nfeval=length(fevalrecline)
      
    }
    if(f<fvalquit){
      if(prtlevel>0)   info=2
      #      return(list(x,f,d,H,iter,info,X,G,w,fevalrec,xrec,Hrec))
    }
    else if(norm(X)>xnormquit){
      if(prtlevel>0) warning("bfgs: norm exceeds specified limit")
      info=3
      #      return(list(x,f,d,H,iter,info,X,G,w,fevalrec,xrec,Hrec))
    }
    if(fail==1){
      if(!quitLSfail){
        if(prtlevel>1) warning("bfgs: continue alhough line search failed")
      }else{
        if(prtlevel>0) warning("bfgs: quit at iteration...")
        info=7
        #        return(list(x,f,d,H,iter,info,X,G,w,fevalrec,xrec,Hrec))
      }
    }
    if(fail==-1){
      if(prtlevel>0) warning("bfgs: f may be unbounded below...")
      info=8
      #      return(list(x,f,d,H,iter,info,X,G,w,fevalrec,xrec,Hrec))
    }
    if(dnorm<=normtol){
      if(prtlevel>0){
        if(nG==1) warning("bfgs: gradient norm below tolerance, quite iteration")
        else warning("bfgs: norm of smallest vector in convex hull of gradients below tolerance, quit...")
      }
      info=0
      #      return(list(x,f,d,H,iter,info,X,G,w,fevalrec,xrec,Hrec))
    }
    #   if(cputime>cpufinish){
    #     if(prtlevel>0) print("bfgs: cpu time limit exceeded, quit iteration")
    #   
    #     info=4
    #     return(list(x,f,d,H,iter,info,X,G,w,fevalrec,xrec,Hrec))
    #   }
    s=alpha*p
    y=g-gprev
    sty=t(s)%*%y
    if(nvec==0){
      if(sty>0){
        if(iter==1 & scale) H=as.numeric(sty/(t(y)%*%y))*H
        rho=as.numeric(1/sty)
        rhoHyst=rho*(H%*%y)%*%t(s)
        H=H-t(rhoHyst)-rhoHyst+rho*s%*%(t(y)%*%rhoHyst)+rho*s%*%t(s)
      }else{
        if(prtlevel>1) warning("bfgs: sty<=0, skipping bfgs update at iteration...")
      }
    }else{
      s=alpha*p
      y=g-gprev
      if(iter<=nvec){
        S=cbind(S,s)
        Y=cbind(Y,y)
      }
      else {
        S=cbind(S[,2:nvec],s)
        Y=cbind(Y[,2:nvec],y)
      }
      if(scale) H=((t(s)%*%y)/(t(y)%*%y))%*%H0
    }
    
  }
  
  return(list(x=x,f=f,d=d,H=H,iter=iter,info=info,X=X,G=G,w=w,fevalrec=fevalrec,xrec=xrec,Hrec=Hrec))
}



















