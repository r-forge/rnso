bfgs=function(pars,options){
  #set options
  x0=options$x0
  x0=t(as.matrix(x0)) #??
  nstart=length(x0)
  fvalquit=options$fvalquit
  xnormquit=options$xnormquit
  prtlevel=options$prtlevel
  #options
  #declarations
  x=list()
  f=c()
  d=list()
  iter=c()
  HH=list()
  H=list()
  info=c()
  X=list()
  G=list()
  w=list()
  fevalrec=c()
  xrec=list()
  Hrec=list()
  
  for(run in 1:nstart){
    tmp=bfgs1run(x0[,run],pars,options)
    #print(tmp)
    x[run]=tmp[[1]]
    f[run]=tmp[[2]]
    d[run]=tmp[[3]]
    HH=tmp[[4]]
    iter[run]=tmp[[5]]
    info[run]=tmp[[6]]
    X[run]=tmp[[7]]
    G[run]=tmp[[8]]
    w[run]=tmp[[9]]
    fevalrec[run]=tmp[[10]]
    xrec[run]=tmp[[11]]
    Hrec[run]=tmp[[12]]
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