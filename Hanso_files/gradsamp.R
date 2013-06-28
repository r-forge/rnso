gradsamp <- function(fn,gr,nvar,x0,f0  =  fn(x0), g0  =  gr(x0), samprad  =  c(1e-3,1e-4,1e-4),maxit  =  100,normtol  =  1e-4, ngrad  =  2, fvalquit  =  -Inf, prtlevel  =  1){
  
  nstart <- ncol(x0)
  f  <-  c()
  x  <-  matrix(NA,nvar,nstart)
  g  <-  matrix(NA,nvar,nstart)
  dnorm  <-  c()
  X  <-  list()
  G  <-  list()
  w  <-  matrix(NA,nvar,nstart)
  for(run in 1:nstart){
    f0 <- fn(x0[,run])
    g0 <- gr(x0[,run])
    
    if((is.nan(f0) | f0  ==  Inf | maxit  == 0)&(prtlevel>0)){
      warning('Gradsamp: function is null or infinite or maxit is zero at initial point')
      f[run]  <-  f0
      x[,run]  <-  x0[,run]
      g[,run]  <-  g0
      dnorm[run]  <-  norm(g0)
      X[[run]]  <-  x[,run]
      G[[run]]  <-  g0
      w[,run]  <-  1
      
    } 
    else {
	tmp  <-  gradsamp1run(fn,gr,x0[,run],f0,g0,samprad,maxit,normtol,ngrad,fvalquit,prtlevel)
	x[,run]  <-  tmp$x
	f[run]  <-  tmp$f
	g[,run]  <-  tmp$g
	dnorm[run]  <-  tmp$dnorm
	X[[run]]  <-  tmp$X
	G[[run]]  <-  tmp$G
	w[,run]  <-  tmp$w
    }
  }
return(list(x  <-  x,f  <-  f,g  <-  g,dnorm  <-  dnorm, X  <-  X, G  <-   G, w  <-  w))
}