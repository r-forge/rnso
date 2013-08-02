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
hanso(fn, gr=NULL, x0 = NULL,dir="forward", nvar = 0, nstart = 10,
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
\item{dir}{direction of gradient approximation to be used. Choose from ("forward","backward","central")}
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
point found by BFGS as its starting point.\cr
\cr
BFGS phase: BFGS is run from multiple starting points, taken from
       the columns of x0, if provided, and otherwise 10 points
       generated randomly. If the termination test was satisfied at the
       best point found by BFGS, HANSO terminates;
       otherwise, it continues to:\cr
       \cr
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
\code{\link{bfgs}}, \code{\link{gradsamp}}, \code{\link{nlcg}},  \code{\link{shor}}
}

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

(res=hanso(fr,grr,nvar=2))

#Using Nesterov's function
#example not run
#fnNesterov1 <- function(x) {
#  n <- length(x)
#  x2 <- x[2:n]; x1 <- x[1:(n-1)]
#  1/4*(x[1]-1)^2 + sum(abs(x2-2*x1^2+1))
#}

#res=hanso(fnNesterov1,nvar=2)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ optimize, nonsmooth }

