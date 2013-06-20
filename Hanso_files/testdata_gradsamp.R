#testdata gradsamp
#rm(list=ls())
#setwd("~/Dropbox/work/gsoc/nso/rnso/pkg")
options=c()
pars=c()
pars$nvar=1
pars$ffgname='fgtest'
options$maxit=1000
options$normtol=1e-6
options$fvalquit=-Inf
options$xnormquit=Inf
options$cpumax=Inf
options$prtlevel=1
options$samprad=c(1e-3,1e-4,1e-5)
options$strongwolfe=0
options$wolfe1=1e-4
options$wolfe2=0.5
options$quitLSfail=1
options$nvec=0
options$H0=1
options$scale=1
options$ngrad=2
options$evaldist=1e-4

x0=c(0.5377  ,  1.8339 ,  -2.2588  ,  0.8622 ,   0.3188  , -1.3077 ,  -0.4336  ,  0.3426 ,   3.5784  ,  2.7694)
options$x0=x0
fgtest=function(x,par=c()){
  f=(x-1)*(x-1)
  g=2*(x-1)
  return(data.frame(f=f,g=g))
}
tmp=fgtest(x0[1])
f0=tmp$f
g0=tmp$g

run=1
choice=1
#gradsampfixed(x0[run], f0, g0, samprad[choice], pars, options)

source('getbundle.R')
source('isnaninf.R')
source('qpspecial.R')
source('gradsampfixed.R')
#getbundle(x0[run], g0, samprad[choice], options$ngrad, pars);