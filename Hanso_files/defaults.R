#default options
bfgs_default  <- function(nvar)
{
  options=c()
  options$nvar=nvar
  options$nstart=10
  options$maxit=1000
  options$normtol=1e-4
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
  options$H0=diag(2)
  options$scale=1
  options$ngrad=2
  options$evaldist=1e-4
  return(options)
}

# source('bfgs.R')
# source('bfgs1run.R')
# source('linesch_ww.R')
# source('hgprod.R')
# source('qpspecial.R')
