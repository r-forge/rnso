f_easy <- function(x){
  fv <- sum(x*x)
  fv <- fv*(1 + 0.1*sin(10*(x[1]+x[2])))
  ifail <- 0
  icount <- 1
  list(fv = fv, ifail = ifail, icount = icount)
}