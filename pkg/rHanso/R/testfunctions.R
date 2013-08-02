##
##  t e s t f u n c t i o n s . R  Non-smooth Test Functions
##


#-- Shor's piecewise quadratic function --------------------------------------
# starting value c(1,1,1,1,1) in [0, 2]^5
# xmin = (1.1243510, 0.9794616, 1.4777077, 0.9202335, 1.1242916)
# fmin = 22.6001621

shor_f <- function(x) {
    stopifnot(is.numeric(x), length(x) == 5)
    x <- as.matrix(c(x))
    A <- matrix(
    c(0,  2,  1,  1,  3,  0,  1,  1,  0,  1,
      0,  1,  2,  4,  2,  2,  1,  0,  0,  1,
      0,  1,  1,  1,  1,  1,  1,  1,  2,  2,
      0,  1,  1,  2,  0,  0,  1,  2,  1,  0,
      0,  3,  2,  2,  1,  1,  1,  1,  0,  0),
    5, 10, byrow = TRUE)
    b <- as.matrix(c(1, 5, 10, 2, 4, 3, 1.7, 2.5, 6, 4.5))
    f <- 0
    for (i in 1:10) {
        d <- b[i] * sum((x - A[, i])^2)
        if (d > f) f <- d
    }
    f
}

shor_g <- function(x) {
    stopifnot(is.numeric(x), length(x) == 5)
    x <- as.matrix(c(x))
    A <- matrix(
    c(0,  2,  1,  1,  3,  0,  1,  1,  0,  1,
      0,  1,  2,  4,  2,  2,  1,  0,  0,  1,
      0,  1,  1,  1,  1,  1,  1,  1,  2,  2,
      0,  1,  1,  2,  0,  0,  1,  2,  1,  0,
      0,  3,  2,  2,  1,  1,  1,  1,  0,  0),
    5, 10, byrow = TRUE)
    b <- as.matrix(c(1, 5, 10, 2, 4, 3, 1.7, 2.5, 6, 4.5))
    f <- 0
    for (i in 1:10) {
        d <- b[i] * sum((x - A[, i])^2)
        if (d > f) {
            f <- d
            k <- i
        }
    }
    2 * b[k] * (x - A[, k])
}





