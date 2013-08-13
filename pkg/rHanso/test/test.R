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
    sum <- sum + 100*(x[i-1]^2-x[i])^2+(x[i-1]-1)^2+90*(x[i+1]^2-x[i+2])^2+(x[i+1]-1)^2+(10*x[i]+x[i+2]-2)^2+((x[i]-x[i+2])^2)/10
  }
  sum
}
res=hanso(fchainwood,x0=c(-3,-1,-3,-1,-2,0))
#compare with optim
resoptim = optim(c(-3,-1,-3,-1,-2,0),fchainwood)

#wood in 4 dimentions
#http://seal.web.cern.ch/seal/documents/minuit/mntutorial.pdf
wood4 <- function(x){
  100*(x[2]-x[1]^2)^2+(x[1]-1)^2+90*(x[4]-x[3]^2)^2+(1-x[3])^2+10.1*((x[2]-1)^2+(x[4]-1)^2)+19.8*(x[2]-1)*(x[4]-1)
}
x0 <- c(-3,-1,-3,-1)
res = hanso(wood4,x0=x0) #100 iterations might not be enough for gradient sampling, for this function.
#try with 2000 iterations for maxitgs, this might take a while
res = hanso(wood4,x0=x0,maxitgs=2000)


#generalized broyden tridiagonal function prob 5
broyden <- function(x){
  n <- length(x)
  x <- c(0,x,0)
  sum <- 0
  for(i in 2:n+1){
    sum <- sum+(abs((3-2*x[i])*x[i]-x[i-1]-x[i+1]+1))^(7/3)
  }
  sum
}
res <- hanso(broyden,nvar=10)
res <- hanso(broyden,x0=rep(-1,10))

#generalized brown function problem 13

genbrown <- function(x){
  n <- length(x)
  k <-floor(n/2)
  sum <-0
  for(j in 1:k){
    i <- 2*j
    sum <- sum + (x[i-1]^2)^(x[i]^2+1)+(x[i]^2)^(x[i-1]^2+1)
  }
  sum
}

fnsphere <- function(x){
  sum(x*x)
}

grsphere <- function(x){
  2*x
}

#nonsmooth version fo brown function
#http://napsu.karmitsa.fi/publications/largetest.pdf (prob2.7)
nsbrown <- function(x)
{
  n <- length(x)
  sum <- 0
  for (i in 1:(n-1)){
    sum <- sum + abs(x[i])^(x[i+1]^2+1)+abs(x[i+1])^(x[i]^2+1)
  }
  sum
}

#http://seal.web.cern.ch/seal/documents/minuit/mntutorial.pdf
powel4 <- function(x){
  (x[1]+10*x[2])^2+5*(x[3]-x[4])^2+(x[2]-2*x[3])^4+10*(x[1]-x[4])^4
}
grpowel4 <- function(x){
  x1 <- 2*(x[1]+10*x[2])+40*(x[1]-x[4])^3
  x2 <- 20*(x[1]+10*x[2])+4*(x[2]-2*x[3])^3
  x3 <- 10*(x[3]-x[4])+8*(2*x[3]-x[2])^3
  x4 <- 10*(x[4]-x[3])+40*(x[4]-x[1])^3
  c(x1,x2,x3,x4)
}
x0=c(3,-1,0,1)
res=hanso(powel4,grpowel4,x0=x0)
