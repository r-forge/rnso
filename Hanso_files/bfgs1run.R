##
##  b f g s . R
##


bfgs1run <- function(fn, gr, x0, H0 = NULL, maxit = 1000,  fvalquit = -Inf,
                     normtol = 1e-6, xnormquit = Inf, evaldist = 1e-4, 
                     ngrad = 0, scale = 1, strongwolfe = 0, 
                     wolfe1 = 1e-4, wolfe2 = 0.5, quitLSfail = TRUE)
{
    n <- length(x0)
    x0 <- as.matrix(x0)
    iter <-0
    # Control parameter checking
    if (is.null(H0)) H0 <- diag(1, n)
    if (ngrad == 0) ngrad = min(100, 2*n, n+10)

    # Utility functions
    is.nainf <- function(x) any(is.na(x) | is.infinite(x))

    # Initializations, all vectors should be column vectors
    x <- x0;    H <- H0
    f <- fn(x); g <- as.matrix(gr(x))
    d <- g
    G <- g;     X <- x
    nG <- 1;    w <- 1
    dnorm <- norm(g, type = "f")

    # Checking error and end conditions
    if (is.nainf(f) || is.nainf(g)) {
        errno <- 5
        mess <- "Function or its gradient not defined at initial iterate."
        return( list(x = c(x), f = f, d = c(d), H = H, iter = iter,
                     info = errno, message = mess, X = X, G = G, w = w))
    }
    if (f < fvalquit) {
        errno <- 2
        mess <- "Function value found below target at initial iterate."
        return( list(x = c(x), f = f, d = c(d), H = H, iter = iter,
                     info = errno, message = mess, X = X, G = G, w = w))
    }
    if (dnorm < normtol) {
        errno <- 0
        mess <- "Tolerance o gradient satisfied at initial iterate."
        return( list(x = c(x), f = f, d = c(d), H = H, iter = iter,
                     info = errno, message = mess, X = X, G = G, w = w))
    }
    if (norm(x, 'f') > xnormquit) {
        errno <- 3
        mess <- "Value of norm(x) exceeds limit at initial iterate."
        return( list(x = c(x), f = f, d = c(d), H = H, iter = iter,
                     info = errno, message = mess, X = X, G = G, w = w))
    }

    #-- Begin main loop ----------------
    for(iter in 1:maxit) {
        #-- Full BFGS, no limited memory version
        p <- -H %*% g
        gtp <- c(t(g) %*% p)
        if (gtp >= 0 || is.na(gtp)) {
            errno <- 6
            mess <- "Found non-descent direction only, will quit."
            return( list(x = c(x), f = f, d = c(d), H = H, iter = iter,
                         info = errno, message = mess, X = X, G = G, w = w))
        }
        gprev <- g  # for BFGS update

        if(strongwolfe){
	  sls <- linesch_sw(fn, gr, x, d = p, f0 = fn(x), grad0 = gr(x),
                          c1 = wolfe1, c2 = wolfe2, fvalquit, prtlevel)
		    alpha <- sls$alpha
		    x <- sls$x
		    f <- sls$f
		    g <- sls$grd
		    fail <- sls$fail
        }
        else{
        wls <- linesch_ww(fn, gr, x, d = p, fn0 = fn(x), gr0 = gr(x),
                          c1 = wolfe1, c2 = wolfe2)
        alpha <- wls$alpha
        x <- as.matrix(wls$xalpha)
        f <- wls$falpha
        g <- as.matrix(wls$galpha)
        fail <- wls$fail
	}
        # discard the saved gradients iff the new point x is not sufficiently
        # close to the previous point and replace them by new gradient
        if (alpha * norm(p, 'f') > evaldist) {
            nG <- 1
            G <- g
            X <- x
        # otherwise add new gradient to set of saved gradients, 
        # discarding oldest if already have ngrad saved gradients
        } else if (nG < ngrad){
            nG <- nG + 1
            G <- cbind(g,G)
            X <- cbind(x,X)
        } else {  # nG = ngrad
            G <- cbind(g, G[, 1:(ngrad-1)])
            X <- cbind(x, X[, 1:(ngrad-1)])
        }

        # compute smallest vector in convex hull of qualifying gradients:
        # reduces to norm of latest gradient if ngrad == 1, and the set must
        # always have at least one gradient: could gain efficiency here by
        # updating previous QP solution
        if (nG > 1) {
            qps <- qpspecial(as.matrix(G))
            w <- qps$x
            d <- qps$d
        } else {
            w <- 1
            d <- g
        }
        dnorm <- sqrt(sum(d*d))

        if (f < fvalquit) {
            errno <- 2
            mess <- "Reached target objective, quit during iteration."
            return( list(x = c(x), f = f, d = c(d), H = H, iter = iter,
                    info = errno, message = mess, X = X, G = G, w = w))
        } else if (norm(x) > xnormquit) {
            errno <- 3
            mess <- "The value of norm(x) exceeds limit during iteration."
            return( list(x = c(x), f = f, d = c(d), H = H, iter = iter,
                    info = errno, message = mess, X = X, G = G, w = w))
        }

        # Line search failed (Wolfe conditions not both satisfied)
        if (fail == 1) {
            if (!quitLSfail) {
                warning("BFGS: Search continued although line search failed.")
            } else {
                errno <- 7
                mess <- "Quit at as line search failed during iteration."
                return( list(x = c(x), f = f, d = c(d), H = H, iter = iter,
                        info = errno, message = mess, X = X, G = G, w = w))
            }
        # Function apparently unbounded from below
        } else if (fail == -1) {
            errno <- 8
            mess <- "Quit as f may be unbounded from below."
            return( list(x = c(x), f = f, d = c(d), H = H, iter = iter,
                         info = errno, message = mess, X = X, G = G, w = w))
        }

        if (dnorm <= normtol) {
            if (nG == 1) {
                mess <- "Gradient norm below tolerance, quit iteration."
            }else {
                mess <- "Norm of smallest vector below tolerance, quit."
            }
            errno <- 0
            return( list(x = c(x), f = f, d = c(d), H = H, iter = iter,
                    info = errno, message = mess, X = X, G = G, w = w))
        }

        # Continue with iteration
        s <- alpha * p
        y <- g - gprev
        sty <- c(t(s) %*% y)

        # Perform rank two BFGS update to the inverse Hessian H
        if (sty > 0) {
            if (iter == 1 && scale) {
                H <- c(sty/(t(y)%*%y)) * H
            }
            
            rho <- 1/sty
            rhoHyst <- rho * (H %*% y) %*% t(s)
            H <- H - t(rhoHyst) - rhoHyst + 
                     rho * s %*% (t(y) %*% rhoHyst) + rho * s %*% t(s)
        } else {
            if (printlevel > 0)
                cat("BFGS: sty<=0 during iteration, skipping bfgs update.\n")
        }

    }  # end for

    mess <- "Maximum number of iterations reached; may not converge."
    errno <- 1
    return( list(x = c(x), f = f, d = c(d), H = H, iter = iter,
            info = errno, message = mess, X = X, G = G, w = w))
}
