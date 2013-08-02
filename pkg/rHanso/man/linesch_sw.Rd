\name{linesch_sw}
\alias{linesch_sw}

\title{
strong Wolfe line search
}
\description{
strong Wolfe line search with cubic interpolation. Can be used by 
using strongwolfe=1 in BFGS. For non smooth functions, this is not
recommended for BFGS. Use weak line search instead. recommended for use 
with CG, where strong Wolfe condition needed for
  convergence analysis
}
\usage{
linesch_sw(fn, gr, x0,dir="forward", d, f0 = fn(x0), 
          grad0 = gr(x0), c1 = 0, c2 = 0.5, fvalquit = -Inf, prtlevel = 0)
}

\arguments{
  \item{fn}{
A function to be minimized. fn(x) takes input as a vector of parameters over which minimization is to take place. fn() returns a scaler.
}
  \item{gr}{
A function to return the gradient for fn(x).
}
  \item{x0}{
initial point
}
\item{dir}{direction of gradient approximation to be used. Choose from ("forward","backward","central")}
  \item{d}{
search direction
}

  \item{f0}{
fn(x0)
}
  \item{grad0}{
gr(x0)
}
  \item{c1}{
Wolfe parameter for the sufficient decrease condition 
}
  \item{c2}{
c2: Wolfe parameter for the WEAK condition on directional derivative
}

  \item{fvalquit}{
quit if f gets below this value.
}
  \item{prtlevel}{
prints messages if this is 1
}
}

\value{
returns a list containing the following fields: \cr
\item{alpha}{steplength satisfying Wolfe conditions}
\item{x}{x0 + alpha*d}
\item{f}{f(x0 + alpha d)}
\item{grad}{(grad f)(x0 + alpha d)}
\item{fail}{0 if both Wolfe conditions satisfied, or falpha < fvalquit
            1 if one or both Wolfe conditions not satisfied but an
               interval was found bracketing a point where both satisfied
           -1 if no such interval was found, function may be unbounded below}
\item{nsteps}{number of steps taken in lszoom}           
}

\references{
Numerical Optimization by 
Jorge Nocedal and   Stephen J. Wright
}
\author{
  Copyright (c) 2010 Michael Overton for Matlab code and documentation,
  with permission converted to R by Abhirup Mallik (and Hans W Borchers).
}
\seealso{
\code{\link{linesch_ww}}
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

res=linesch_sw(fr,grr,c(-1.2,1),c(1,1),dir="f")
}
\keyword{ optimize }
