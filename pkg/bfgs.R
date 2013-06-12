bfgs=function(pars,options){
  #set options
  x0=options$x0
  nstart=length(x0)
  fvalquit=options$fvalquit
  xnormquit=options$xnormquit
  prtlevel=options$prtlevel
  #options
  for(run in 1:nstart){
    tmp=bfgs1run(x0[,run],pars,options)
    x[,run]=tmp$x
    f[run]=tmp$f
    d[,run]=tmp$d
    HH=tmp$H
    iter[run]=tmp$iter
    info[run]=tmp$info
    X[run]=tmp$X
    G[run]=tmp$G
    w[run]=tmp$w
  }
  H[run]=(HH+t(HH))/2
  #cputime break condition
  if(nstart==1){
    H=H[1]
    fevalrec=fevalrec[1]
    xrec=xrec[1]
    Hrec=Hrec[1]
    X=X[1]
    G=G[1]
    w=w[1]
  }
  return(list(x,f,d,H,iter,X,G,fevalrec,xrec,Hrec))
}