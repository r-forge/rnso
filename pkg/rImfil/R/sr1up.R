sr1up <- function(x, xc, sgrad, gc, hess, alist) {
  n <- length(x)
  pr <- diag(alist)
  y <- sgrad - gc
  s <- x - xc
  z <- y - hess%*%s
  
  y <- pr%*%y
  z <- pr%*%z
  
  if (sum(z*s) != 0) {
    ptst <- sum(z*(hess%*%z))+(sum(z*z)^2)/sum(z*s)
    if (ptst > 0) {
      hess <- pr %*% hess %*% pr + (z %*% t(z))/sum(z*s)
      hess <- diag(n) - pr + hess
    }
  }
  if (kappa(hess) > 1e6) hess <- diag(n)
  hess
}