grad_nso <- function(fn, x0, dir = c("forward", "backward", "central"), heps = .Machine$double.eps^(1/2), 
    ...) {
    if (!is.numeric(x0)) 
        stop("Argument 'x0' must be a numeric value.")
    
    fun <- match.fun(fn)
    fn <- function(x) fun(x, ...)
    if (length(fn(x0)) != 1) 
        stop("Function 'f' must be a univariate function of n variables.")
    
    dir <- match.arg(dir)
    
    n <- length(x0)
    hh <- rep(0, n)
    gr <- matrix(NA, nrow = , ncol = 1)
    for (i in 1:n) {
        hh[i] <- heps
        if (dir == "forward") {
            gr[i] <- (fn(x0 + hh) - fn(x0))/heps
        } else if (dir == "backward") {
            gr[i] <- (fn(x0) - fn(x0 - hh))/heps
        } else if (dir == "central") {
            gr[i] <- (fn(x0 + hh) - fn(x0 - hh))/(2 * heps)
        } else {
            stop("Direction must be one of 'forward', 'backward', or 'central'.")
        }
        hh[i] <- 0
    }
    return(gr)
} 
