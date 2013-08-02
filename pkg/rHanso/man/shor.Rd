\name{shor}
\alias{shor}
\title{
Shor's R Algorithm
}
\description{
Shor's R algorithm for unconstrained minimization of smooth and non smooth functions.
}
\usage{
shor(fn, gr,dir="forward", nvar = 0, nstart = 10, 
    x0 = NULL, maxit = 1000, fvalquit = -Inf,
    beta = 0.5, normtol = 1e-06, xnormquit = Inf,
    evaldist = 1e-04, ngrad = 0, rescale = 0,
    strongwolfe = 0, useprevstep = 0, wolfe1 = 1e-04,
    wolfe2 = 0.5, quitLSfail = TRUE, prtlevel = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fn}{
A function to be minimized. fn(x) takes input as a vector of parameters over which minimization is to take place. fn() returns a scaler.
}
  \item{gr}{
Optional input. A function to return the gradient for fn(x). If not provided, will be calculated internally.
}
\item{dir}{direction of gradient approximation to be used. Choose from ("forward","backward","central")}
  \item{nvar}{
Number of parameters that fn() takes.
}
  \item{nstart}{
Number of initial guesses. Default is 10.
}
  \item{x0}{
Matrix, with dimension (nvar x nstart), each column represent an initial guess.
}
  \item{maxit}{
maximum number of iterations. 
}

  \item{fvalquit}{
quit if f drops below this value.
}
  \item{beta}{
the key parameter, between 0 and 1, or nan to minimize
          the error in the secant equation each iteration (default 1/2)
          beta = 1: steepest descent;  beta -> 0: conjugate gradient
           (beta = 0 will give divide by 0)
            note: beta = 1 - gamma, where gamma is notation in IMAJNA paper
}
  \item{normtol}{
termination tolerance on d: smallest vector in
           convex hull of up to ngrad gradients
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
  \item{rescale}{
1 if rescale B to have inf norm 1 every iteration
}
  \item{strongwolfe}{
if this is 1, strong line search is used, otherwise, weak line search is used. Default is 0.
}
  \item{useprevstep}{
if 1, line search is initialized with previous steps. 1 seemed to perform better, but hard to justify this rationally.
}
  \item{wolfe1}{
wolfe line search parameter 1.
}
  \item{wolfe2}{
wolfe line search parameter 2.
}
  \item{quitLSfail}{
1 if want to quit when line search fails. 0 otherwise.
}  \item{prtlevel}{
if 1, prints in code messages.
}
}

\value{
Returns a list containing the following fields:
\item{x}{a matrix with k'th column containing final value of x obtained from k'th column of x0.}
\item{f}{a vector of final obtained minimum values of fn() at the initial points.}
\item{g}{each column is the gradient at the corresponding column of x}
\item{B}{list of the final shor matrices}
\item{frec}{record of function evaluations}
\item{betarec}{record of beta}
\item{xrec}{record of iterates}
\item{svrec}{record of singular value's of the Hessian}
}

\author{
  Copyright (c) 2010 Michael Overton for Matlab code and documentation,
  with permission converted to R by Abhirup Mallik (and Hans W Borchers).
  }
\seealso{
\code{\link{hanso}}, \code{\link{gradsamp}}, \code{\link{nlcg}},  \code{\link{bfgs}}
}
\examples{
fnNesterov2 <- function(x) {
 n <- length(x)
 x2 <- x[2:n]; x1 <- x[1:(n-1)]
 1/4*abs((x[1]-1)) + sum(abs(x2-2*abs(x1)+1))
}
res=shor(fnNesterov2,nvar=5)
res$f
}
\keyword{ optimize }
