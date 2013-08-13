nlcg <- function(fn, gr = NULL, nvar = 0, nstart = 10, x0 = NULL, upper = 1, lower = 0, H0 = NULL, maxit = 1000, 
    fvalquit = -Inf, prtlevel = 0, version = "C", normtol = 1e-06, xnormquit = Inf, evaldist = 1e-04, ngrad = 0, 
    scale = 1, wolfe1 = 1e-04, wolfe2 = 0.5, quitLSfail = TRUE, strongwolfe = 1) {
    if (!is.null(x0)) {
        
        if (class(x0) == "numeric") {
            x0 <- t(x0)
            nstart <- 1
            nvar = length(x0)
        } else if (class(x0) == "matrix") {
            nvar <- nrow(x0)
            nstart <- ncol(x0)
        } else stop("unknown initial value matrix, please enter a numeric vector or matrix")
    } else {
        nstart <- 10
        M <- matrix(runif(nvar * nstart), nrow = nvar, ncol = nstart)
        x0 <- (upper - lower) * M + lower
    }
    
    
    x <- matrix(NA, nvar, nstart)
    f <- c()
    g <- matrix(NA, nvar, nstart)
    frec <- list()
    alpharec <- list()
    message <- list()
    for (run in 1:nstart) {
        tmp <- nlcg1run(fn, gr, x0[, run], H0, maxit, fvalquit, strongwolfe, version, prtlevel, normtol, 
            xnormquit, evaldist, ngrad, scale, wolfe1, wolfe2, quitLSfail)
        # print(tmp)
        x[, run] <- tmp$x
        f[run] <- tmp$f
        g[, run] <- tmp$g
        frec[[run]] <- tmp$frec
        alpharec[[run]] <- tmp$alpharec
        message[[run]] <- tmp$message
    }
    list(x = x, f = f, g = g, frec = frec, alpharec = alpharec, message = message)
} 
