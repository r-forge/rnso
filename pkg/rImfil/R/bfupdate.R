bfupdate <- function(x, xc, sgrad, gc, hess, alist) {
  n <- length(x)
  pr <- diag(alist)
  y <- sgrad - gc
  s <- x - xc
  x <- hess %*% s
  y <- as.numeric(pr) * y
  if (sum(y*s) > 0) {
    hess <- as.numeric(pr) * hess * as.numeric(pr) + (y%*%t(y)/sum(t*s)) - as.numeric(pr)*(z%*%t(z)/sum(s*z))*as.numeric(pr)
    hess <- diag(n) - pr + hess
  }
  if (kappa(hess) > 1e6){
    hess <- diag(n)
  }
  hess
}