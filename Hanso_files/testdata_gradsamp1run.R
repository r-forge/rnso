source('gradsampfixed.R')
source('qpspecial.R')
source('isnaninf.R')
source('linesch_ww.R')
source('getbundle.R')
source('gradsamp1run.R')
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
x0 <- c(1.2,3.8)
tmp <- gradsamp1run(fr,grr,x0)
gradsampfixed(fr,grr,x0,samprad=1e-3)