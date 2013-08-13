shor1run <- function(fn, gr, x0, maxit = 1000, fvalquit = -Inf, beta = 0.5, normtol = 1e-06, xnormquit = Inf, 
    evaldist = 1e-04, ngrad = 0, rescale = 0, strongwolfe = 0, useprevstep = 0, wolfe1 = 1e-04, wolfe2 = 0.5, 
    quitLSfail = TRUE, prtlevel = 1) {
    
    is.nainf <- function(x) any(is.na(x) | is.infinite(x))
    
    nvar <- length(x0)
    x <- x0
    f <- fn(x)
    g <- gr(x)
    
    
    alpha <- 1
    gnorm <- sqrt(sum(g * g))
    n <- nvar
    B <- diag(n)
    I <- diag(n)
    frec <- c()
    fevalrecall <- list()
    betarec <- c()
    xrec <- matrix(NaN, n, maxit)
    svrec <- matrix(NaN, n, maxit)
    if (is.nainf(f)) {
        mess <- "shor: f is infinite or na at initial iteration."
        if (prtlevel) 
            cat(mess)
        return(list(x = c(x), f = f, g = g, B = B, frec = frec, fevalrecall = fevalrecall, betarec = betarec, 
            xrec = xrec, svrec = svrec))
    } else if (gnorm < normtol) {
        mess <- "shor: gradient less than the tolerance."
        if (prtlevel) 
            cat(mess)
        return(list(x = c(x), f = f, g = g, B = B, frec = frec, fevalrecall = fevalrecall, betarec = betarec, 
            xrec = xrec, svrec = svrec))
    } else if (f < fvalquit) {
        mess <- "shor: f is below target objective."
        if (prtlevel) 
            cat(mess)
        return(list(x = c(x), f = f, g = g, B = B, frec = frec, fevalrecall = fevalrecall, betarec = betarec, 
            xrec = xrec, svrec = svrec))
    }
    
    for (iter in 1:maxit) {
        
        p <- -B %*% (t(B) %*% g)
        if (useprevstep) 
            p <- alpha * p
        gtp <- sum(g * p)
        
        if (gtp >= 0 || is.na(gtp)) {
            mess <- paste("shor: not descent direction, quit at iter=", iter, "\n")
            if (prtlevel) 
                cat(mess)
            return(list(x = c(x), f = f, g = g, B = B, frec = frec, fevalrecall = fevalrecall, betarec = betarec, 
                xrec = xrec, svrec = svrec))
        }
        gprev <- g
        if (strongwolfe) {
            nfeval <- NaN
            res_sw <- linesch_sw(fn, gr, x, p, f, g, wolfe1, wolfe2, fvalquit, prtlevel)
            alpha <- res_sw$alpha
            x <- res_sw$x
            f <- res_sw$f
            g <- res_sw$grd
            fail <- res_sw$fail
            
            if (wolfe2 == 0) {
                increase <- 1e-04 * alpha + 1e-08
                x <- x + increase * p
                f <- fn(x)
                g <- gr(x)
            }
            fevalrec <- NaN
        } else {
            res_ww <- linesch_ww(fn, gr, x, p, f, g, wolfe1, wolfe2, fvalquit, prtlevel)
            alpha <- res_ww$alpha
            x <- res_ww$xalpha
            f <- res_ww$falpha
            g <- res_ww$galpha
            fail <- res_ww$fail
            fevalrec <- res_ww$fevalrec
        }
        
        gnorm <- sqrt(sum(g * g))
        xrec[, iter] <- x
        fevalrecall[[iter]] <- fevalrec
        svrec[, iter] <- svd(B)$d
        frec[iter] <- f
        
        if (f < fvalquit) {
            mess <- paste("shor: reached target objective, quit at iter=", iter, "\n")
            if (prtlevel) 
                cat(mess)
            return(list(x = c(x), f = f, g = g, B = B, frec = frec, fevalrecall = fevalrecall, betarec = betarec, 
                xrec = xrec, svrec = svrec))
        }
        if (fail == 1) {
            mess <- paste("shor: linesearch failed, at iter=", iter, "\n")
            if (prtlevel) 
                cat(mess)
            if (quitLSfail) {
                return(list(x = c(x), f = f, g = g, B = B, frec = frec, fevalrecall = fevalrecall, betarec = betarec, 
                  xrec = xrec, svrec = svrec))
            }
        } else if (fail == -1) {
            mess <- paste("shor: function could be unbounded below, quit at iter=", iter, "\n")
            if (prtlevel) 
                cat(mess)
            return(list(x = c(x), f = f, g = g, B = B, frec = frec, fevalrecall = fevalrecall, betarec = betarec, 
                xrec = xrec, svrec = svrec))
        }
        if (gnorm <= normtol) {
            mess <- paste("shor: gradient norm below tolerance, quit at iter=", iter, "\n")
            if (prtlevel) 
                cat(mess)
            return(list(x = c(x), f = f, g = g, B = B, frec = frec, fevalrecall = fevalrecall, betarec = betarec, 
                xrec = xrec, svrec = svrec))
        }
        y <- g - gprev
        r <- t(B) %*% y
        xi <- r/sqrt(sum(r * r))
        if (is.na(beta)) {
            y <- g - gprev
            s <- alpha * p
            Hy <- B %*% r
            HHy <- B %*% (t(B) %*% Hy)
            betasq <- sum(s * Hy)/(sum(y * HHy))
            
            if (betasq > 1) 
                betasq <- 1 else if (betasq < 0) 
                betasq <- 0 else if (is.na(betasq)) 
                betasq <- 0.25
            
            beta <- sqrt(betasq)
        }
        betarec[iter] <- beta
        B <- B %*% (I + (beta - 1) * (xi %*% t(xi)))
        if (rescale) 
            B <- B/norm(B, "I")
    }
    mess <- paste("shor: maximum iterations reached, quitting", "\n")
    if (prtlevel) 
        cat(mess)
    return(list(x = c(x), f = f, g = g, B = B, frec = frec, fevalrecall = fevalrecall, betarec = betarec, 
        xrec = xrec, svrec = svrec))
}






















 
