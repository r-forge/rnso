f=1.2148
g=-2.2044
x =
-0.1022
p=2.2044
pars=c()
pars$nvar=1
pars$fgname='fgtest'
wolfe1=1e-4
wolfe2=0.500
fvalquit=-Inf
prtlevel=1

fgtest=function(x,par=c()){
  f=(x-1)*(x-1)
  g=2*(x-1)
  return(data.frame(f=f,g=g))
}
source('linesch_ww.R')
#linesch_ww(x, f, g, p, pars, wolfe1, wolfe2, fvalquit, prtlevel)

#testing linesch_ww coded by Hans:


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

linesch_ww(fr,grr,c(-1.2,1),c(1,1))
