rm(list=ls())
source('bfgs.R')
source('bfgs1run.R')
source('linesearch.R')
source('hgprod.R')
source('qpspecial.R')
#source('defaults.R')
library(numDeriv)
source('linesch_sw.R')
source('lszoom.R')
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
res=bfgs(fr,grr,nvar=2)
res=bfgs(fr,grr,nvar=2,strongwolfe=1)

##nesterov's function

fnNesterov1 <- function(x) {
  n <- length(x)
  x2 <- x[2:n]; x1 <- x[1:(n-1)]
  1/4*(x[1]-1)^2 + sum(abs(x2-2*x1^2+1))
}

grNest <-function(x){
  grad(fnNesterov1,x)}

(res=bfgs(fnNesterov1,grNest,nvar=5))

## just any max(abs) function
testfn <- function(x) {
  return ( max(abs(x[1]),norm(as.matrix(x[-1],'F'))))
}

grtest <-function(x){
  grad(testfn,x)}
#change nvar as you like
(res1=bfgs(testfn,grtest,nvar=5))
