getbundle=function(x,g,samprad,N,pars){
  m=length(x)
  #declare empty matrices
  xbundle[,1]=x
  gbundle[,1]=g
  for(k in 2:N){
    xpert=x+samprad%*%(runif(m)-0.5)
    tmp=fgtest(xpert,pars)
    f=tmp$f
    grad=tmp$g
    count=0
    while(isnaninf(f) | isnaninf(grad)){
      xpert=(x+xpert)/2
      f=fgtest(xpert,pars)
      count=count+1
    }
    xbundle[,k]=xpert
    gbundle[,k]=xpert
  }
  return(list(xbundle,gbundle))
}