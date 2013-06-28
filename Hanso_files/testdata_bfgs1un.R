# options=c()
# pars=c()
# pars$nvar=2
# options$maxit=1000
# options$normtol=1e-4
# options$fvalquit=-Inf
# options$xnormquit=Inf
# options$cpumax=Inf
# options$prtlevel=1
# options$samprad=c(1e-3,1e-4,1e-5)
# options$strongwolfe=0
# options$wolfe1=1e-4
# options$wolfe2=0.5
# options$quitLSfail=1
# options$nvec=0
# options$H0=diag(2)
# options$scale=1
# options$ngrad=2
# options$evaldist=1e-4
# 
# options$x0=rbind(c( 1.3703 ,  -0.1022  ,  0.3192 ,  -0.8649 ,  -0.1649 ,   1.0933 ,  -0.8637  , -1.2141,  -0.0068  , -0.7697),c( -1.7115 ,  -0.2414 ,   0.3129 ,  -0.0301  ,  0.6277  ,  1.1093  ,  0.0774  , -1.1135,  1.5326   , 0.3714))
# # fgtest=function(x,par=c()){
# #   f=(x-1)*(x-1)
# #   g=2*(x-1)
# #   return(data.frame(f=f,g=g))
# # }

fr <- function(x) {   ## Rosenbrock Banana function
  x1 <- x[1]
  x2 <- x[2]
  100 * (x2 - x1 * x1)^2 + (1 - x1)^2
}
grr <- function(x) { ## Gradient of 'fr'
  x1 <- x[1]
  x2 <- x[2]
  c(-400 * x1 * (x2 - x1 * x1) - 2 * (1 - x1),
    200 *      (x2 - x1 * x1))
}
source('linesch_ww.R')
source('qpspecial.R')
source('isnaninf.R')
source('bfgs1run.R')
res <- bfgs1run(fr,grr,2,c(-1.2,1))

