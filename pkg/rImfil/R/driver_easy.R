driver_easy <- function(){
  bounds <- matrix(c(-1,1,-1,1),2,2)
  budget <- 40
  tmp <- imfil(x0, f_easy, budget, bounds)
  x <- tmp$x
  histout <- tmp$histout
  list(x = x, histout = histout)
}