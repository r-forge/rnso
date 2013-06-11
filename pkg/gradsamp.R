gradsamp=function(pars,options){
  options=setdefaults(pars,options)
  options=setx0(pars,options)
  x0=options$x0
  nstart=length(x0)
  nvar=pars$nvar
  maxit=options$maxit
  #cpufinish=cputime+options$cpumax
  ptm=proc.time()
  prtlevel=options$prtlevel
  if(is.na(options$samprad)) options$samprad=c(1e-4,1e-5,1e-6)
  if(is.na(options$ngrad)) options$ngrad=min(100,2*nvar,nvar+100)
  for(run=1:nstart){
    #if(prtlevel>0 & nstart >1) print("gradsamp: starting point ")
    tmp=fgtest(x0[,run],pars)
    f0=tmp$f
    g0=tmp$g
#     f[run]=f0
#     x[,run]=x0[,run]
#     dnorm[,run]=x[,run]
#     X=list() ##?
#     G=list() ##?
#     w[,run]=1##
     tmp=gradsamp1run(x0[,run],f0,g0,pars,options) ## catch returns
    if((proc.time()-ptm)[3]>1000000) break
  }
  ## nstart==1
}