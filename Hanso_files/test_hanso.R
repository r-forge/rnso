source('hanso.R')
source('postprocess.R')
source('gradsampfixed.R')
source('qpspecial.R')
source('isnaninf.R')
source('linesch_ww.R')
source('getbundle.R')
source('gradsamp1run.R')
source('gradsamp.R')
source('bfgs.R')
source('bfgs1run.R')
source('hgprod.R')
source('qpspecial.R')
library(numDeriv)

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

res=hanso(fr,grr,2)

fnNesterov1 <- function(x) {
  n <- length(x)
  x2 <- x[2:n]; x1 <- x[1:(n-1)]
  1/4*(x[1]-1)^2 + sum(abs(x2-2*x1^2+1))
}

grNest <-function(x){
  grad(fnNesterov1,x)}

(res=hanso(fnNesterov1,grNest,nvar=2))

## just any max(abs) function
testfn <- function(x) {
  return ( max(abs(x[1]),norm(as.matrix(x[-1],'F'))))
}

grtest <-function(x){
  grad(testfn,x)}
#change nvar as you like
(res1=hanso(testfn,grtest,nvar=5))

fg <- function(x)
{
  
  x1 <- x[1]; x2 <- x[2]
  res<- 100*(x2 - x1*x1)^2 + (1-x1)^2
  return(res)
}

gr <- function(x) {
  x1 <- x[1]; x2 <- x[2]
  c(-400*x1*(x2 - x1*x1)-2*(1-x1), 200*(x2 - x1*x1))
}

tmp=hanso(fg,gr,nvar=2)