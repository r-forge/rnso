\name{imfil}
\alias{imfil}
\title{
Minimization of noisy functions subject to explicit bound
 constraints and hidden constraints
}
\description{
This is an implementation of Kelly's Implicit filtering algorithm for minimizing functions under contrained and bounded space.
}
\usage{
imfil(x0, fn, budget, bounds, options = imfil_optset())
}

\arguments{
  \item{x0}{
Initial iterate
}
  \item{fn}{
Objective function
}
  \item{budget}{
maximum cost
}
  \item{bounds}{
The interval c(low, high) for x.
}
  \item{options}{
options for the function. The defaults are set by imfil_optset()
}
}

\details{
Implicit filtering solves bounded constrained optimization problems, where the goal is to minimize the objective function f subject to the condition that x remains in the feasible region. Implicit filtering is a sampling method. The optimization is controlled only by evaluating f at a cluster of points in the sample space. That evaluation determines the next cluster.
}
\value{
\item{x}{estimated minimizer.}
\item{histout}{iteration history. Updated after each nonlinear iteration. 
}
\item{complete_history}{complete evaluation history}
}
\references{
Iterative Methods for Optimization, by C. T. Kelly
}
\author{
Converted by Abhirup Mallik and Hans Borchers from Matlab Code by C.T. Kelly
}

\examples{
\dontrun{
fn <- f_easy
bounds <- matrix(c(-1,-1,1,1),2,2)
budget <- 40
x0 <- matrix(c(0.5,0.5))
results <- imfil(x0, fn, budget, bounds, options = imfil_optset())
}
}
\keyword{ optimize }

