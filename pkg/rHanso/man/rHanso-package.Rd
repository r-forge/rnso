\name{rHanso-package}
\alias{rHanso-package}
\alias{rHanso}
\docType{package}
\title{
Hybrid Algorithm for Non-Smooth Optimization (HANSO)
}
\description{
This package implements HANSO algorithm for optimization of smooth and nonsmooth functions. It uses
both BFGS and gradient sampling methods. It also provides the users the option to use BFGS, gradient sampling and wolfe line search (both strong and weak).
}
\details{
\tabular{ll}{
Package: \tab rnso\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-07-17\cr
License: \tab MIT\cr
}
The main functions of this package are:\cr
hanso(): Hanso algorithm\cr
bfgs(): Minimiza function using BFGS\cr
gradsamp(): Gradient sampling algorithm for non smooth (and smooth) functions\cr
linesch_sw(): Strong wolfe line search\cr
linesch_ww(): Weak wolfe line search
}
\author{
Abhirup Mallik, University of Minnesota, Hans Werner Borchers

Maintainer: Abhirup Mallik <malli066@umn.edu>
}
\references{
A.S. Lewis and M.L. Overton, Nonsmooth Optimization via BFGS, 2008. \cr
    J.V. Burke, A.S. Lewis and M.L. Overton, A Robust Gradient Sampling 
     Algorithm for Nonsmooth, Nonconvex Optimization
      SIAM J. Optimization 15 (2005), pp. 751-779
}

\keyword{ package }
\examples{
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
#BFGS
(res=bfgs(fr,grr,nvar=2)) # using all the default parameters, using weak wolfe search
#HANSO
(res=hanso(fr,grr,nvar=2))
#Gradient sampling
x0 <- matrix(rnorm(12),2,6)
(tmp <- gradsamp(fr,grr,2,x0))
#Strong line search
(res=linesch_sw(fr,grr,c(-1.2,1),c(1,1)))
#Weak line search
linesch_ww(fr,grr,c(-1.2,1),c(1,1))

}
