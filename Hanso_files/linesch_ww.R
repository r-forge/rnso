linesch_ww=function(x0,f0,grad0,d,pars,c1,c2,fvalquit,prtlevel){
  ## check for valid input ##
  fgname=pars$fgname
  alpha=0
  xalpha=x0
  falpha=f0
  gradalpha=grad0
  beta=Inf
  gradbeta=NA*rep(1,length(x0))
  g0=t(as.matrix(grad0))%*%as.matrix(d)
  if(g0>=0) print("linesch_ww_mod:Warning, not a descent direction")
  dnorm=norm(as.matrix(d),type="1")
  if(dnorm==0) print("d is zero")
  t=1
  nfeval=0
  nbisect=0
  nexpand=0
  nbisectmax=max(30,round(log2(1e5*dnorm)))
  nexpandmax=max(30,round(log2(1e5/dnorm)))
  done=0
  fevalrec=c()
  while (!done){
    x=x0+t*d
    nfeval=nfeval+1
    tmp=fgtest(x,par)
    f=tmp$f
    grad=tmp$g
    fevalrec[nfeval]=f
    if(f<fvalquit){
      fail=0
      alpha=t
      xalpha=x
      falpha=f
      gradalpha=grad
      return(list(alpha, xalpha, falpha, gradalpha, fail, beta, gradbeta, fevalrec))
    }
    gtd=t(grad)%*%as.matrix(d)
    if(f>= f0+c1*t*g0 | is.na(f)){
      beta=t
      gradbeta=grad
    }
    else if(gtd<=c2*g0 | is.na(gtd))
    {
      alpha=t
      xalpha=x
      falpha=f
      gradalpha=grad
    }else{
      fail=0
      alpha=t
      xalpha=x
      falpha=f
      gradalpha=grad
      beta=t
      gradbeta=grad
      return(list(alpha, xalpha, falpha, gradalpha, fail, beta, gradbeta, fevalrec))
    }
    if(beta<Inf){
      if(nbisect<nbisectmax){
	nbisect=nbisect+1
	t=(alpha+beta)/2
      }
      else done=1
      
    }else{
      if(nexpand<nexpandmax){
	nexpand=nexpand+1
	t=2*alpha
      }
      else done=1
    }
  }
  if(beta==Inf){
    fail=-1
    if(prtlevel>1){
      print("Line search failed to bracket")
      print("wolfe conditions, function may be unbounded below")
    }}
    else{
      fail=1
      if(prtlevel>1)
	{ print("Line search failed to satisfy weak wolfe conditions")
	  print("although point satisfying conditions was bracketed")
	}
  }
}























