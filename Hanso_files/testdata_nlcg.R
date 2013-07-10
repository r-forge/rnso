source('lszoom.R')
source('linesch_sw.R')
source('linesearch.R')
source('nlcg1run.R')
source('nlcg.R')
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
res=nlcg(fr,grr,2)
x0=matrix(c(0.2,1.3,0.9,.24),2,2)
res=nlcg(fr,grr,2,nstart=2,x0=x0,strongwolfe=0)
res=nlcg(fr,grr,2,nstart=2,x0=x0)