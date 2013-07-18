linesch_ww <-
function( fn, gr=NULL, x0, d, fn0 = fn(x0), gr0 = gr(x0),
                        c1 = 0, c2 = 0.5,
                        fvalquit = -Inf, prtlevel = 0 ) {
    
  if(is.null(gr)){
    gr <- function(x){
      grad_nso(fn,x,dir="forward")
    }
  }
    stopifnot(is.numeric(x0), is.numeric(d))
    if (c1 < 0 || c1 > c2 || c2 > 1)    # 0 <= c1 <= c2 <= 1/n
        stop("Arguments 'c1','c2' must satisfy: 0 <= c1 <= c2 <= 1/n.")
    n <- length(x0)

    #-- steplength parameters
    alpha  <- 0         # lower bound on steplength conditions
    xalpha <- x0
    falpha <- fn0
    galpha <- gr0       # need to pass grad0, not grad0'*d, in case line search fails
    beta   <- Inf       # upper bound on steplength satisfying weak Wolfe conditions
    gbeta  <- rep(NA, n)

    g0 <- sum(gr0 * d)
    if (g0 >= 0)
        if (prtlevel > 0)
            mess <- paste("Linesearch: Argument 'd' is not a descent direction.")
    dnorm <- sqrt(sum(d * d))
    if (dnorm == 0)
        stop("Linesearch: Argument 'd' must have length greater zero.")

    t <- 1              # important to try step length one first
    nfeval   <- 0
    nbisect  <- 0
    nexpand  <- 0
    nbisectmax <- max(30, round(log2(1e5*dnorm)))   # allows more if ||d|| big
    nexpandmax <- max(10, round(log2(1e5/dnorm)))   # allows more if ||d|| small

    #-- main loop ----------------------
    fevalrec <- c()
    fail <- 0
    done <- FALSE
    while (!done) {
        x <- x0 + t*d
        fun <- fn(x)
        grd <- gr(x)
        nfeval <- nfeval + 1
        fevalrec <- c(fevalrec, fun)
        if (fun < fvalquit) {
            return( list(alpha = t, xalpha = x, falpha = fun, galpha = grd,
                        fail = fail, beta = beta, gbeta = gbeta, fevalrec = fevalrec))
        }

        gtd <- sum(grd * d)
        if (fun >= fn0 + c1*t*g0 || is.na(fun)) {# first condition violated
            beta  <- t
            gbeta <- grd
        } else if (gtd <= c2*g0 || is.na(gtd)) {# second condition violated   
            alpha  <- t
            xalpha <- x
            falpha <- fun
            galpha <- grd
        } else {# both conditions satisfied           
            return( list(alpha = t, xalpha = x, falpha = fun, galpha = grd,
                         fail = fail, beta = t, gbeta = grd, fevalrec = fevalrec))
        }

        # set up next function evaluation
        if (beta < Inf) {
            if (nbisect < nbisectmax) {
                nbisect <- nbisect + 1
                t <- (alpha + beta)/2           # bisection
            } else {
                done <- TRUE
            }
        } else {
            if (nexpand < nexpandmax) {
                nexpand <- nexpand + 1
                t <- 2*alpha                    # still in expansion mode
            } else {
                done <- TRUE
            }
        }
    } # end while

    # Wolfe conditions not satisfied; there are two cases:
    if (is.infinite(beta)) {# minimizer never bracketed
        fail <- -1
        if (prtlevel > 0)
            mess <- paste("Linesearch: Function may be unbounded from below.")
    } else {# point satisfying Wolfe conditions bracketed
        fail <- 1
        if (prtlevel > 0)
            mess <- paste("Linesearch: Failed to satisfy weak Wolfe conditions.")
    }
    
    return( list(alpha = t, xalpha = x, falpha = fun, galpha = grd,
                 fail = fail, beta = t, gbeta = grd, fevalrec = fevalrec) )
}
