\name{bfgs}
\alias{bfgs}
\title{Optimization using BFGS
}
\description{
BFGS is normally used for optimizing smooth, not necessarily convex, 
  functions, for which the convergence rate is generically superlinear.
   But it also works very well for functions that are nonsmooth at their  
   minimizers, typically with a linear convergence rate and a final 
   inverse Hessian approximation that is very ill conditioned, as long 
   as a weak Wolfe line search is used. This version of BFGS will work
   well both for smooth and nonsmooth functions and has a stopping 
   criterion that applies for both cases, described above.
}
\usage{
bfgs(fn, gr, nvar=0, nstart = 10, x0 = NULL, upper = 1,
      lower = 0,  maxit = 1000, normtol = 1e-06,
      fvalquit = -Inf, xnormquit = Inf, nvec = 0,
      prtlevel = 1, strongwolfe = 0, wolfe1 = 1e-04,
      wolfe2 = 0.5, quitLSfail = 1, ngrad = 0, evaldist = 1e-04,
      H0 = NULL, scale = 1)
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
  \item{nstart}{
Number of initial guesses. Default is 10.
}
  \item{x0}{
Matrix, with dimension (nvar x nstart), each column represent an initial guess.
}
\item{upper}{upper bound for the initial value}
\item{lower}{lower bound for the initial value}
  \item{maxit}{
maximum number of iterations. 
}
  \item{normtol}{
termination tolerance on d: smallest vector in
           convex hull of up to ngrad gradients
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
}

\value{
Returns a list containing the following fields:
\item{x}{a matrix with k'th column containing final value of x obtained from k'th column of x0.}
\item{f}{a vector of final obtained minimum values of fn() at the initial points.}
\item{d}{the final smallest vector in the convex hull of the saved gradient.}
\item{H}{final BFGS inverse hessian approximation}
\item{iter}{number of iterations}
\item{message}{any warnings / messages for checking the termination condition.}
\item{X}{iterates, where saved gradients were evaluated.}
\item{G}{saved gradients used for computation.}
\item{w}{weights giving the smallest vector.}

}
\references{
A.S. Lewis and M.L. Overton, Nonsmooth Optimization via BFGS, 2008.
}
\author{
  Copyright (c) 2010 Michael Overton for Matlab code and documentation,
  with permission converted to R by Abhirup Mallik (and Hans W Borchers).
}
\seealso{
\code{\link{hanso}}, \code{\link{gradsamp}}, \code{\link{nlcg}},  \code{\link{shor}}
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
res=bfgs(fr,grr,nvar=2) # using all the default parameters, using weak wolfe search
res
res=bfgs(fr,grr,nvar=2,strongwolfe=1) # using strong wolfe search
res
# Nesterov's function in 5 dimention
# example not run
#fnNesterov1 <- function(x) {
#  n <- length(x)
#  x2 <- x[2:n]; x1 <- x[1:(n-1)]
#  1/4*(x[1]-1)^2 + sum(abs(x2-2*x1^2+1))
#}

#res=bfgs(fnNesterov1,nvar=5)
#res

}

\keyword{optimization}

