nsf <- function(x) {
  + f1 <- x[1]^2 + x[2]^2
  + f2 <- x[1]^2 + x[2]^2 + 10 * (-4*x[1] - x[2] + 4)
  + f3 <- x[1]^2 + x[2]^2 + 10 * (-x[1] - 2*x[2] + 6)
  + max(f1, f2, f3)
  + }

fn=function(x){
  x*x*(2+sin(pi/x))
}

fr <- function(x) {  
  x1 <- x[1]
  x2 <- x[2]
  (x2 - x1 * x1)^2 
}
grr <- function(x) { 
  x1 <- x[1]
  x2 <- x[2]
  c(-4 * x1 * (x2 - x1 * x1),
    2 *      (x2 - x1 * x1))
}


fn=function(x){
  s=0
  for (i in 1:1){
    s=s+(x[i+1]-2*x[i]*x[i]+1)
    }
  (x[1]-1)*(x[1]-1)/4+s
  }
