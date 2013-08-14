\name{nlcg}
\alias{nlcg}
\title{
Nonlinear Conjugate Gradient minimization
}
\description{
   NLCG (Nonlinear Conjugate gradient)
   is intended for minimizing smooth, not necessarily convex, functions.
}
\usage{
nlcg(fn, gr = NULL,nvar = 0, nstart = 10,
      x0 = NULL,upper = 1, lower = 0,  H0 = NULL, maxit = 1000, fvalquit = -Inf,
      prtlevel = 0, version = "C", normtol = 1e-06,
      xnormquit = Inf, evaldist = 1e-04, ngrad = 0, scale = 1,
      wolfe1 = 1e-04, wolfe2 = 0.5, quitLSfail = TRUE, strongwolfe = 1)
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
  \item{H0}{
Initial inverse hessian approximation. Default is NULL}
  \item{scale}{1 to scale H0 at first iteration, 0 otherwise.
}
  \item{maxit}{
maximum number of iterations. 
}
  \item{fvalquit}{
quit if f drops below this value.
}
  \item{prtlevel}{
prints messages if this is 1
}
  \item{version}{
'P' for Polak-Ribiere-Polyak (not recommended: fails on hard problems)
           'F' for Fletcher-Reeves (not recommended: often stagnates)
           'C' for Polak-Ribiere-Polyak Constrained by Fletcher-Reeves 
              (recommended, combines advantages of 'P' and 'F'; default)
           'S' for Hestenes-Stiefel (not recommended)
           'Y' for Dai-Yuan (allows weak Wolfe line search, see below)
           'Z' for Hager-Zhang
           '-' for Steepest Descent (for comparison)
}
  \item{normtol}{
termination tolerance on d: smallest vector in
           convex hull of up to options.ngrad gradients
          (default: 1e-6) 
}
  \item{xnormquit}{
quit if norm(x) drops below this.
}
  \item{evaldist}{
the gradients used in the termination test
           qualify only if they are evaluated at points  approximately 
           within distance options.evaldist of x
}
  \item{ngrad}{
number of gradients willing to save and use in
           solving QP to check optimality tolerance on smallest vector in
           their convex hull; see also next two options 
}
  \item{wolfe1}{
wolfe line search parameter 1.
}
  \item{wolfe2}{
wolfe line search parameter 2.
}
  \item{quitLSfail}{
if 1, quits when line search is failed. If 0, continues.
}
  \item{strongwolfe}{
if 1, strong wolfe line search is used.
}
}
\details{
   NLCG (Nonlinear Conjugate gradient)
   is intended for minimizing smooth, not necessarily convex, functions.
   The Fletcher-Reeves version (version='F') is globally convergent in
   theory but often stagnates in practice.  The Polak-Ribiere version
   (version='P') works better in practice but its search direction
   may not even be a descent direction and it may not converge.  
   The 'C' version combines the best of both.  It is Polak-Ribiere 
   constrained by Fletcher-Reeves, typically behaving as well as or better
   than PR in practice but with the same global convergence guarantee as FR.
   The Hestenes-Stiefel version (version='S') also often fails.
   The Dai-Yuan (version='Y') is the first to allow use of a weak Wolfe
   line search.  The Hager-Zhang (version='Z') is the newest and is also
   promising.

}
\value{
Returns a list containing the following fields:
\item{x}{a matrix with k'th column containing final value of x obtained from k'th column of x0.}
\item{f}{a vector of final obtained minimum values of fn() at the initial points.}
\item{g}{each column is the gradient at the corresponding column of x}
\item{frec}{recorded f values for all iterates}
\item{alpharec}{record of steps taken in line search}
}
\references{
 Reference: second edition of Nocedal and Wright, Chapter 5, plus papers
   by Dai-Yuan and Hager-Zhang.
}
\author{
  Copyright (c) 2010 Michael Overton for Matlab code and documentation,
  with permission converted to R by Abhirup Mallik (and Hans W Borchers).
}

\seealso{
\code{\link{hanso}}, \code{\link{gradsamp}}, \code{\link{bfgs}},  \code{\link{shor}}
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
res=nlcg(fr,grr,nvar=2)
## Following examples are not run, please uncomment them to run.
#x0=matrix(c(0.2,1.3,0.9,.24),2,2)
#res=nlcg(fr,grr,nvar=2,nstart=2,x0=x0,strongwolfe=0) #using weak wolfe
}
\keyword{ optimize }
