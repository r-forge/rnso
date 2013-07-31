fnNesterov2 <- function(x) {
 n <- length(x)
 x2 <- x[2:n]; x1 <- x[1:(n-1)]
 1/4*abs((x[1]-1)) + sum(abs(x2-2*abs(x1)+1))
}

fnNesterov1 <- function(x) {
  n <- length(x)
  x2 <- x[2:n]; x1 <- x[1:(n-1)]
  1/4*(x[1]-1)^2 + sum(abs(x2-2*x1^2+1))
}


#grNest1 <-function(x){
#  grad(fnNesterov1,x)}

#grNest2 <-function(x){
#  grad(fnNesterov2,x)}

fr <- function(x) {   ## Rosenbrock Banana function
  x1 <- x[1]
  x2 <- x[2]
  100 * (x2 - x1 * x1)^2 + (1 - x1)^2
}
grr <- function(x) { ## Gradient of 'fr'
  x1 <- x[1]
  x2 <- x[2]
  c(-400 * x1 * (x2 - x1 * x1) - 2 * (1 - x1),
    200 *      (x2 - x1 * x1))
}

(res=hanso(fr,grr,nvar=2))

#chained rosenbrock function prob 2
frosen <- function(x){
  n=length(x)
  x1=x[-1]
  x=x[-n]
  sum(100*(x^2-x1^2)^2+(x-1)^2)
}

res=hanso(frosen,nvar=6)

#chained wood function prob 2
fchainwood <- function(x){
  n <- length(x)
  sum <-0
  for(j in 1:((n-2)/2)){
    i <- 2*j
    sum <- sum + 100*(x[i-1]^2-x[i])^2+(x[i-1]-1)^2+90*(x[i+1]^2-x[i+2])^2
  }
  sum
}
res=hanso(fchainwood,nvar=6)


