\name{hanso}
\alias{hanso}
\title{
HANSO: Hybrid Algorithm for Nonsmooth Optimization
}
\description{
Minimization algorithm intended for nonsmooth, nonconvex functions,
   but also applicable to functions that are smooth, convex or both.

}
\usage{
hanso(fn, gr = NULL, x0 = NULL,upper = 1, lower = 0, nvar = 0, nstart = 10,
      maxit = 1000, maxitgs=100, normtol = 1e-06, fvalquit = -Inf,
      xnormquit = Inf, nvec = 0, prtlevel = 1, strongwolfe = 0,
      wolfe1 = 1e-04, wolfe2 = 0.5, quitLSfail = 1,
      ngrad = min(100, 2 * nvar, nvar + 10), evaldist = 1e-04,
      H0 = diag(nvar), scale = 1, samprad = c(1e-04, 1e-05, 1e-06))
}
\arguments{
  \item{fn}{
A function to be minimized. fn(x) takes input as a vector of parameters over which minimization is to take place. fn() returns a scaler.
}
  \item{gr}{
Optional input. A function to return the gradient for fn(x). If not provided, will be calculated internally.
}
  \item{x0}{
Matrix, with dimension (nvar x nstart), each column represent an initial guess.
}
 \item{upper}{upper bound for the initial values}
 \item{lower}{lower bound for initial values.}
  \item{nvar}{
Number of parameters that fn() takes.
}
  \item{nstart}{
Number of initial guesses. Default is 10.
}
  \item{maxit}{
maximum number of iterations for bfgs. 
}
  \item{maxitgs}{
  maximum number of iterations for gradient sampling. Warning!: Gradient Sampling is expensive.
  }
  \item{normtol}{
termination tolerance on d: smallest vector in
           convex hull of up to options.ngrad gradients
          (default: 1e-6) 
}
  \item{fvalquit}{
quit if f drops below this value.
}
  \item{xnormquit}{
quit if norm(x) drops below this.
}
  \item{nvec}{
0 for saving bfgs values.
}
  \item{prtlevel}{
prints messages if this is 1
}
  \item{strongwolfe}{
if this is 1, strong line search is used, otherwise, weak line search is used. Default is 0.
}
  \item{wolfe1}{
wolfe line search parameter 1.
}
  \item{wolfe2}{
wolfe line search parameter 2.
}
  \item{quitLSfail}{
1 if want to quit when line search fails. 0 otherwise.
}
  \item{ngrad}{
number of gradients willing to save and use in
           solving QP to check optimality tolerance on smallest vector in
           their convex hull; see also next two options 
}
  \item{evaldist}{
the gradients used in the termination test
           qualify only if they are evaluated at points  approximately 
           within distance options.evaldist of x
}
  \item{H0}{
Initial inverse hessian approximation. Default is NULL}
  \item{scale}{1 to scale H0 at first iteration, 0 otherwise.
}
  \item{samprad}{
vector of sampling radii. For gradient sampling.
}
}
\details{
It is a two phase process, BFGS phase works for smooth functions, gradient
sampling phase works for non smooth ones. Gradient sampling uses the minimum 
point found by BFGS as its starting point.

BFGS phase: BFGS is run from multiple starting points, taken from
       the columns of x0, if provided, and otherwise 10 points
       generated randomly. If the termination test was satisfied at the
       best point found by BFGS, HANSO terminates;
       otherwise, it continues to:

      Gradient sampling phases: 3 gradient sampling phases are run from 
       lowest point found, using sampling radii: 
       10*evaldist, evaldist, evaldist/10
}
\value{
Returns a list containing the following fields:
\item{x}{a matrix with k'th column containing final value of x obtained from k'th column of x0.}
\item{f}{a vector of final obtained minimum values of fn() at the initial points.}
\item{loc}{local optimality certificate, list with 2 fields:
      dnorm: norm of a vector in the convex hull of gradients of the 
          function evaluated at and near x 
      evaldist: specifies max distance from x at which these gradients 
       were evaluated. The smaller loc$dnorm and loc$evaldist are, the
       more likely it is that x is an approximate local minimizer.}
\item{H}{final BFGS inverse hessian approximation}
\item{X}{iterates, where saved gradients were evaluated.}
\item{G}{saved gradients used for computation.}
\item{w}{weights giving the smallest vector.}
       
}
\references{
A.S. Lewis and M.L. Overton, Nonsmooth Optimization via BFGS, 2008. \cr
    J.V. Burke, A.S. Lewis and M.L. Overton, A Robust Gradient Sampling 
     Algorithm for Nonsmooth, Nonconvex Optimization
      SIAM J. Optimization 15 (2005), pp. 751-779
}
\author{
  Copyright (c) 2010 Michael Overton for Matlab code and documentation,
  with permission converted to R by Abhirup Mallik (and Hans W Borchers).
  \email{malli066@umn.edu}
}
\seealso{
\code{\link{shor}}, \code{\link{gradsamp}}
}

\examples{
## Nesterov function in 2 dimensions
nesterov_f <- function(x) {
    x1 <- x[1]; x2 <- x[2]
    1/4*(x1-1)^2 + abs(x2 - 2*x1^2+1)
 }

nesterov_g <- function(x) {     # analytical gradient
    g <- matrix(NA, 2, 1)
    x1 <- x[1]; x2 <- x[2]
    sgn <- sign(x2 - 2*x1^2 + 1)
    if (sgn != 0) {
        g[1] <- 0.5*(x1-1) - sgn*4*x1
        g[2] <- sgn
    }
    g
 }

 hanso(nesterov_f, nesterov_g, nvar = 2)
 
 # Hanso: Best value found by BFGS =  0.1401677 
 # gradsamp: not descent direction, quit at iter= 93 
 # gradsamp: not descent direction, quit at iter= 9 
 # gradsamp: not descent direction, quit at iter= 1 
 # Best value found by Gradient Sampling =  1.065261e-09 
 # $x
 #           [,1]
 # [1,] 0.9999712
 # [2,] 0.9998849
 #
 # $f
 # [1] 1.065261e-09
 # 
 # ...

\dontrun{
## Shor's piecewise quadratic function
hanso(shor_f, shor_g, nvar=5)

# Hanso: Best value found by BFGS =  24.99301 
# gradsamp: not descent direction, quit at iter= 100 
# gradsamp: not descent direction, quit at iter= 38 
# gradsamp: not descent direction, quit at iter= 7 
# Best value found by Gradient Sampling =  22.60016 
# $x
#           [,1]
# [1,] 1.1243379
# [2,] 0.9794677
# [3,] 1.4777124
# [4,] 0.9202101
# [5,] 1.1242998
# 
# $f
# [1] 22.60016
# 
# ...
 }
}

\keyword{ optimize }
\keyword{ nonsmooth }
