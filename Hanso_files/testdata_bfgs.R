rm(list=ls())
source('bfgs.R')
source('bfgs1run.R')
source('linesch_ww.R')
source('hgprod.R')
source('qpspecial.R')
source('defaults.R')
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
bfgs(fr,grr,2)
