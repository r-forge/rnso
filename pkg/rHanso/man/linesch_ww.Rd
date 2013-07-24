\name{linesch_ww}
\alias{linesch_ww}
\title{
Weak wolfe line search
}
\description{
Line search enforcing weak Wolfe conditions, suitable
           for minimizing both smooth and nonsmooth functions. Intended to be called from Hanso, but can be used by itself also.
}
\usage{
linesch_ww(fn, gr, x0, d, fn0 = fn(x0), gr0 = gr(x0), c1 = 0, c2 = 0.5, fvalquit = -Inf, prtlevel = 0)
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
  \item{d}{
search direction
}
  \item{fn0}{
fn(x0)
}
  \item{gr0}{
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
\details{

 The weak Wolfe line search is far less complicated that the standard 
 strong Wolfe line search that is discussed in many texts. It appears
 to have no disadvantages compared to strong Wolfe when used with
 Newton or BFGS methods on smooth functions, and it is essential for 
 the application of BFGS or bundle to nonsmooth functions as done in HANSO.
 However, it is NOT recommended for use with conjugate gradient methods,
 which require a strong Wolfe line search for convergence guarantees.
 Weak Wolfe requires two conditions to be satisfied: sufficient decrease
 in the objective, and sufficient increase in the directional derivative
 (not reduction in its absolute value, as required by strong Wolfe).

 There are some subtleties for nonsmooth functions.  In the typical case
 that the directional derivative changes sign somewhere along d, it is
 no problem to satisfy the 2nd condition, but descent may not be possible
 if the change of sign takes place even when the step is tiny. In this
 case it is important to return the gradient corresponding to the positive 
 directional derivative even though descent was not obtained. On the other 
 hand, for some nonsmooth functions the function decrease is steady
 along the line until at some point it jumps to infinity, because an
 implicit constraint is violated.  In this case, the first condition is
 satisfied but the second is not. All cases are covered by returning
 the end points of an interval [alpha, beta] and returning the function 
 value at alpha, but the gradients at both alpha and beta. 

 The assertion that [alpha,beta] brackets a point satisfying the
 weak Wolfe conditions depends on an assumption that the function 
 f(x + td) is a continuous and piecewise continuously differentiable 
 function of t, and that in the unlikely event that f is evaluated at
 a point of discontinuity of the derivative, g'*d, where g is the 
 computed gradient, is either the left or right derivative at the point
 of discontinuity, or something in between these two values.

 For functions that are known to be nonsmooth, setting the second Wolfe
 parameter to zero makes sense, especially for a bundle method, and for
 the Shor R-algorithm, for which it is essential.  However, it's not
 a good idea for BFGS, as for smooth functions this may prevent superlinear 
 convergence, and it can even make trouble for BFGS on, e.g., 
 f(x) = x_1^2 + eps |x_2|, when eps is small.
}

\value{
returns a list containing:
\item{alpha}{steplength satisfying weak Wolfe conditions if one was found,
             otherwise left end point of interval bracketing such a point
             (possibly 0)}
\item{xalpha}{x0 + alpha*d}
\item{falpha}{f(x0 + alpha d)}
\item{gradalpha}{(grad f)(x0 + alpha d)}
\item{fail}{0 if both Wolfe conditions satisfied, or falpha < fvalquit
            1 if one or both Wolfe conditions not satisfied but an
               interval was found bracketing a point where both satisfied
           -1 if no such interval was found, function may be unbounded below}
\item{beta}{same as alpha if it satisfies weak Wolfe conditions,
             otherwise right end point of interval bracketing such a point
             (inf if no such finite interval found)}
\item{gradbeta}{(grad f)(x0 + beta d) (this is important for bundle methods)
             (vector of nans if beta is inf)}
\item{fevalrec}{record of function evaluations}
}

\author{
  Copyright (c) 2010 Michael Overton for Matlab code and documentation,
  with permission converted to R by Abhirup Mallik (and Hans W Borchers).
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

linesch_ww(fr,grr,c(-1.2,1),c(1,1))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ optimization }

