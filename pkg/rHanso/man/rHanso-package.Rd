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
Version: \tab 0.2-1\cr
Date: \tab 2013-07-26\cr
License: \tab GPL (>= 3)\cr
}
The main functions of this package are:
\itemize{
\item  \code{\link{hanso}}: Hanso algorithm
\item \code{\link{bfgs}}: Minimiza function using BFGS
\item \code{\link{gradsamp}}: Gradient sampling algorithm for non smooth (and smooth) functions
\item \code{\link{linesch_sw}}: Strong wolfe line search
\item \code{\link{linesch_ww}}: Weak wolfe line search
\item \code{\link{nlcg}}: Non Linear Conjugate gradient minimization
\item \code{\link{shor}}: An inplementation of Shor's R algorithm
}
}
\author{
Abhirup Mallik and Hans W Borchers
Maintainer: Abhirup Mallik \email{malli066@umn.edu}
}
\references{
A.S. Lewis and M.L. Overton, Nonsmooth Optimization via BFGS, 2008.\cr
URL: http://www.cs.nyu.edu/overton/papers/

J.V. Burke, A.S. Lewis and M.L. Overton, A Robust Gradient Sampling 
Algorithm for Nonsmooth, Nonconvex Optimization
SIAM J. Optimization 15 (2005), pp. 751-779
}

\keyword{ package }
\examples{
# Typical examples will be provided later.
}
