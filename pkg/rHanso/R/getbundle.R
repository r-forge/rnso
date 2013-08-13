getbundle <- function(fn, gr, x, g = gr(x), samprad, N) {
    is.naninf <- function(v) {
        any(is.na(v)) || any(is.infinite(v)) || any(is.nan(v))
    }
    m <- length(x)
    # declare empty matrices
    xbundle <- matrix(NA, m, N)
    gbundle <- matrix(NA, m, N)
    xbundle[, 1] <- x
    gbundle[, 1] <- g
    for (k in 2:N) {
        xpert <- x + samprad * (runif(m) - 0.5)  #samprad is a scaler here
        f <- fn(xpert)
        grd <- gr(xpert)
        count <- 0
        while (is.nainf(f) || is.nainf(grd)) {
            xpert <- (x + xpert)/2
            f <- fn(xpert)
            grd <- gr(xpert)
            count <- count + 1
            if (count > 100) 
                stop("getbundle: too many contractions needed to find finite f and g")
        }
        xbundle[, k] <- xpert
        gbundle[, k] <- grd
    }
    return(list(xbundle = xbundle, gbundle = gbundle))
} 
