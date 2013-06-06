linesch_ww=function(x0,f0,grad0,d,pars,c1,c2,fvalquit,prtlevel){
  ## check for valid input ##
  fgname=pars$fgname
  alpha=0
  xalpha=x0
  falpha=f0
  gradalpha=grad0
  beta=inf
  gradbeta=nan*rep(1,length(x0))
  g0=t(grad0)%*%d
  if(g0>=0) print("error1")
  dnorm=norm(d,type="1")
  if(dnorm==0) print("error2")
  t=1
  nfeval=0
  nbisect=0
  nexpand=0
  nbisectmax=max(30,round(log2(1e5*dnorm)))
  nexpandmax=max(30,round(log2(1e5/dnorm)))
  done=0
  while (?){
    x=x0+t*d
    nfeval=nfeval+1
    tmp=eval(fgname,x,pars)
    f=tmp$f
    grad=tmp$g
    ??
    if(f<fvalquit){
      fail=0
      alpha=t
      xalpha=x
      falpha=f
      gradalpha=grad
      return(c(alpha, xalpha, falpha, gradalpha, fail, beta, gradbeta, fevalrec))
    }
    gtd=t(grad)%*%d
    if(f>= f0+c1*t*g0 | is.nan(f)){
      beta=t
      gradbeta=grad
    }
    if(gtd<=c2*g0*|is.nan(gtd))
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
      return(c(alpha, xalpha, falpha, gradalpha, fail, beta, gradbeta, fevalrec))
    }
    if(beta<inf){
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
  if(beta==inf){
    fail=-1
    if(prtlevel>1){
      print("Line search failed to bracket")
      print("wolfe con...")
    }
    else{
    fail=1
    if(prtlevel>1)
    {
      print("Line search failed to bracket")
      print("wolfe con...")
    }
  }
}























