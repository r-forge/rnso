\name{testfunctions}
\alias{shor_f}
\alias{shor_g}
\title{
  Non-smooth Test Functions
}
\description{
  Examples of test functions for non-smooth optimization.
}
\usage{
shor_f(x)
shor_g(x)
}
\arguments{
  \item{x}{numerical vector of arguments for test functions.}
}
\details{
  The following non-smooth test functions are provided:

  \bold{Shor's function}\cr
  Shor's piecewise quadratic function \code{shor_f} of five variables with 
  gradient function \code{shor_g}.
  In \code{[0, 2]^5} it has its minimum \code{fmin = 22.6001621} at 
  \code{xmin = c(1.1243510, 0.9794616, 1.4777077, 0.9202335, 1.1242916)}.
}
\value{
  All test functions will return a numerical scalar value. The gradient 
  functions will return a column vector of the same length as \code{x}.
}
\examples{
    x5 <- runif(5)
    shor_f(x5)
    shor_g(x5)
}

