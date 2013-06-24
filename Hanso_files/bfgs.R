bfgs=function(fn,gr,nvar,options=bfgs_default(nvar)){
  #set options
  options$x0 <- matrix(rnorm(options$nvar*options$nstart),options$nvar,options$nstart)
  nstart=options$nstart
  fvalquit=options$fvalquit
  xnormquit=options$xnormquit
  prtlevel=options$prtlevel
  #options
  #declarations
  x=list()
  f=c()
  d=list()
  iter=c()
  #HH=list()
  H=list()
  info=c()
  X=list()
  G=list()
  w=list()
  fevalrec=c()
  xrec=list()
  Hrec=list()
  
  for(run in 1:nstart){
    tmp=bfgs1run(fn,gr,options$x0[,run],options)
    #print(tmp)
    x[[run]]=tmp$x
    f[run]=tmp$f
    d[[run]]=tmp$d
    HH=as.matrix(tmp$H)
    iter[run]=tmp$iter
    info[run]=tmp$info
    X[[run]]=tmp$X
    G[[run]]=tmp$G
    w[[run]]=tmp$w
    fevalrec[run]=tmp$fevalrec
    xrec[[run]]=tmp$xrec
    Hrec[[run]]=tmp$Hrec
  }
  H[[run]]=(HH+t(HH))/2
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