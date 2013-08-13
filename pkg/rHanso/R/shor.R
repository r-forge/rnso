shor <- function(fn, gr = NULL, nvar = 0, nstart = 10, x0 = NULL, upper = 1, lower = 0, maxit = 1000, fvalquit = -Inf, 
    beta = 0.5, normtol = 1e-06, xnormquit = Inf, evaldist = 1e-04, ngrad = 0, rescale = 0, strongwolfe = 0, 
    useprevstep = 0, wolfe1 = 1e-04, wolfe2 = 0.5, quitLSfail = TRUE, prtlevel = 1) {
    
    if (!is.null(x0)) {
        
        if (class(x0) == "numeric") {
            x0 <- matrix(x0)
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
    B <- list()
    frec <- list()
    fevalrec <- list()
    betarec <- list()
    xrec <- list()
    svrec <- list()
    
    for (run in 1:nstart) {
        res <- shor1run(fn, gr, x0[, run], maxit, fvalquit, beta, normtol, xnormquit, evaldist, ngrad, rescale, 
            strongwolfe, useprevstep, wolfe1, wolfe2, quitLSfail, prtlevel)
        x[, run] <- res$x
        f[run] <- res$f
        g[, run] <- res$g
        B[[run]] <- res$B
        frec[[run]] <- res$frec
        fevalrec[[run]] <- res$fevalrecall
        betarec[[run]] <- res$betarec
        xrec[[run]] <- res$xrec
        svrec[[run]] <- res$svrec
        
    }
    
    return(list(x = x, f = f, g = g, B = B, frec = frec, fevalrec = fevalrec, betarec = betarec, xrec = xrec, 
        svrec = svrec))
    
} 
