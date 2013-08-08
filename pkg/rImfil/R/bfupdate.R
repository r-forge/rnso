bfupdate <- function(x, xc, sgrad, gc, hess, alist) {
  n <- length(x)
  pr <- diag(alist)
  y <- sgrad - gc
  s <- x - xc
  x <- hess %*% s
  y <- pr %*% y
  if (sum(y*s) > 0) {
    hess <- pr %*% hess %*% pr + (y%*%t(y)/sum(t*s)) - pr%*%(z%*%t(z)/sum(s*z))%*%pr
    hess <- diag(n) - pr + hess
  }
  if (kappa(hess) > 1e6){
    hess <- diag(n)
  }
  hess
}