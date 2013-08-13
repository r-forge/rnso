\name{gradsamp}
\alias{gradsamp}
\title{
Gradient sampling algorithm for Non-Smooth Non-Convex Minimization
}
\description{
Intended for nonconvex functions that are continuous everywhere and for 
   which the gradient can be computed at most points, but which are known 
   to be nondifferentiable at some points, typically including minimizers.
}
\usage{
gradsamp(fn, gr, nvar, x0,upper = 1, lower = 0, 
        f0 = fn(x0), g0 = gr(x0), 
        samprad = c(1e-04, 1e-05, 1e-06), 
        maxit = 1000, normtol = 1e-06, 
        ngrad = min(100, 2 * nvar, nvar + 10),
        fvalquit = -Inf, prtlevel = 1)
}
\arguments{
  \item{fn}{
A function to be minimized. fn(x) takes input as a vector of parameters over which minimization is to take place. fn() returns a scaler.
}
  \item{gr}{
A function to return the gradient for fn(x).
}
  \item{nvar}{
Number of parameters that fn() takes.
}
  \item{x0}{
Initial guess of x.
}
\item{upper}{upper bound for the initial value}
\item{lower}{lower bound for the initial value}
  \item{f0}{
fn(x0). Optional, is taken as input to reduce computation.
}
  \item{g0}{
gr(x0). Optional, is taken as input to reduce computation.
}
  \item{samprad}{
vector of sampling radii. 
}
  \item{maxit}{
number of maximum iterations
}
  \item{normtol}{
stopping tolerance for norm(d)
}
  \item{ngrad}{
number of sampled gradients per iterate.
}
  \item{fvalquit}{
quit, if f reaches this value.
}
  \item{prtlevel}{
  prints messages if this is 1
}
}

\value{
Returns a list containing the following fields:
\item{x}{a matrix with k'th column containing final value of x obtained from k'th column of x0.}
\item{f}{a vector of final obtained minimum values of fn() at the initial points.}
\item{g}{matrix with k'th column as the gradient of fn at x[,k]}
\item{dnorm}{dnorm(k) is the norm of a vector in the convex hull of 
       gradients of the function evaluated at points near x(:,k)}
\item{X}{a list containing points where gradients were evalueated}
\item{G}{a list with gradients of correspoding elements from X}
\item{w}{vectors of weights needed for dnorm}
\item{message}{any warnings / use it to know about termination condition}
}
\references{
J.V. Burke, A.S. Lewis and M.L. Overton,
    A Robust Gradient Sampling Algorithm for Nonsmooth, Nonconvex Optimization
    SIAM J. Optimization, 2005
}
\author{
  Copyright (c) 2010 Michael Overton for Matlab code and documentation,
  with permission converted to R by Abhirup Mallik (and Hans W Borchers).
}
\seealso{
\code{\link{hanso}}, \code{\link{bfgs}}, \code{\link{nlcg}},  \code{\link{shor}}
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
#x0 <- matrix(rnorm(12),2,6)
x0 <- c(1.2,1)
tmp <- gradsamp(fr,grr,x0=x0,maxit=100)
tmp

#using Nesterov's function
# example not run
#fnNesterov1 <- function(x) {
#  n <- length(x)
#  x2 <- x[2:n]; x1 <- x[1:(n-1)]
#  1/4*(x[1]-1)^2 + sum(abs(x2-2*x1^2+1))
#}

#res=gradsamp(fnNesterov1,nvar=2,rep(1,2),maxit=10)
#res
}

\keyword{ optimize, nonsmooth}

