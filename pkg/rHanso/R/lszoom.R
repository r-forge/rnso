lszoom <- function(fn, gr, lo, hi, flo, fhi, glo, ghi, f0, g0, x0, d, c1, c2, prtlevel) {
    
    
    fail <- 0
    lo2 <- hi
    flo2 <- fhi
    glo2 <- ghi
    nsteps <- 0
    
    while ((lo != hi) && (nsteps < 50)) {
        nsteps <- nsteps + 1
        bisect <- (lo + hi)/2
        interp <- cubic_interp(lo, lo2, flo, flo2, glo, glo2)
        if (inside(interp, lo, bisect)) 
            atry <- interp else atry <- bisect
        
        xtry <- x0 + atry * d
        ftry <- fn(xtry)
        gradtry <- gr(xtry)
        gtry <- sum(gradtry * d)
        
        if (ftry > f0 + c1 * atry * g0 || ftry >= flo) {
            hi <- atry
            lo2 <- hi
            flo2 <- ftry
            glo2 <- gtry
        } else {
            if (abs(gtry) <= -c2 * g0) {
                alpha <- atry
                x <- xtry
                f <- ftry
                grd <- gradtry
                return(list(alpha = alpha, x = x, f = f, grd = grd, fail = fail, nsteps = nsteps))
            }
            if (gtry * (hi - lo) >= 0) {
                hi <- lo
            }
            lo2 <- lo
            flo2 <- flo
            glo2 <- glo
            lo <- atry
            flo <- ftry
            glo <- gtry
        }
        
    }
    
    if (prtlevel > 1) 
        cat("linesch_sw: failed to satisfy wolfe conditions, lszoom loop ran for 50 times\n")
    alpha <- atry
    x <- xtry
    f <- ftry
    grd <- gradtry
    fail <- 1
    list(alpha = alpha, x = x, f = f, grd = grd, fail = fail, nsteps = nsteps)
    
}

cubic_interp <- function(x1, x2, f1, f2, g1, g2) {
    eta <- g1 + g2 - 3 * (f1 - f2)/(x1 - x2)
    if (eta^2 > g1 * g2) {
        gamma <- sign(x2 - x1) * sqrt(eta^2 - g1 * g2)
    } else gamma <- NaN
    xmin <- x2 - (x2 - x1) * (g2 + gamma - eta)/(g2 - g1 + 2 * gamma)
    xmin
}

inside <- function(x, a, b) {
    isin <- 0
    if (is.na(x)) 
        return(isin)
    if (a <= b) {
        if (x >= a & x <= b) 
            isin <- 1
    } else {
        if (x >= b & x <= a) 
            isin <- 1
    }
    isin
}
 
