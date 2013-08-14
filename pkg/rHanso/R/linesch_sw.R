linesch_sw <- function(fn, gr, x0, d, f0 = fn(x0), grad0 = gr(x0),
             c1 = 0, c2 = 0.5, fvalquit = -Inf, prtlevel = 0) {
  
  
  stopifnot(is.numeric(x0), is.numeric(d))
  if (c1 < 0 || c1 > c2 || c2 > 1) 
    stop("Arguments 'c1','c2' must satisfy: 0 <= c1 <= c2 <= 1/n.")
  n <- length(x0)
  g0 <- sum(grad0 * d)
  if (g0 >= 0) 
    if (prtlevel > 0) 
      warning("Linesearch: Argument 'd' is not a descent direction.")
  dnorm <- sqrt(sum(d * d))
  if (dnorm == 0) 
    stop("Linesearch: Argument 'd' must have length greater zero.")
  if (fvalquit != -Inf) 
    stop("Linesearch: option 'fvalquit' has not yet been implemented.")
  
  
  old <- 0
  fold <- f0
  gold <- g0
  new <- 1
  
  nexpand <- max(50, round(log2(dnorm)))
  for (k in 1:nexpand) {
    xnew <- x0 + new * d
    fnew <- fn(xnew)
    gradnew <- gr(xnew)
    gnew <- sum(gradnew * d)
    if (fnew > f0 + c1 * new * g0 | ((fnew >= fold) & k > 1)) {
      tmp <- lszoom(fn, gr, old, new, fold, fnew, gold, gnew, f0, g0, x0, d, c1, c2, prtlevel)
      alpha <- tmp$alpha
      x <- tmp$x
      f <- tmp$f
      grd <- tmp$grd
      fail <- tmp$fail
      nsteps <- tmp$nsteps
      return(list(alpha = alpha, x = x, f = f, grd = grd, fail = fail, steps = nsteps))
    }
    
    if (abs(gnew) <= -c2 * g0) {
      alpha <- new
      x <- xnew
      f <- fnew
      grd <- gradnew
      fail <- 0
      nsteps <- k
      return(list(alpha = alpha, x = x, f = f, grd = grd, fail = fail, steps = nsteps))
    }
    if (gnew >= 0) {
      tmp <- lszoom(fn, gr, new, old, fnew, fold, gnew, gold, f0, g0, x0, d, c1, c2, prtlevel)
      alpha <- tmp$alpha
      x <- tmp$x
      f <- tmp$f
      grd <- tmp$grd
      fail <- tmp$fail
      nsteps <- tmp$nsteps
      return(list(alpha = alpha, x = x, f = f, grd = grd, fail = fail, steps = nsteps))
    }
    
    old <- new
    fold <- fnew
    gold <- gnew
    new <- 2 * new
  }
  mess <- "linesch_sw: minimizer was not bracketed, function could be unbounded below"
  alpha <- new
  x <- xnew
  f <- fnew
  grd <- gnew
  fail <- -1
  nsteps <- 0
  
  list(alpha = alpha, x = x, f = f, grd = grd, fail = fail, nsteps = nsteps, message = mess)
  
} 
