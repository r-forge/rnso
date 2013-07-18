gradsamp <-
function(fn,gr,nvar,x0,f0  =  fn(x0), g0  =  gr(x0), samprad  =  c(1e-4,1e-5,1e-6),maxit  =  1000,normtol  =  1e-6, ngrad  =  min(100,2*nvar,nvar+10), fvalquit  =  -Inf, prtlevel  =  1){
  x0=as.matrix(x0)
  nstart <- ncol(x0)
  f  <-  c()
  x  <-  matrix(NA,nvar,nstart)
  g  <-  matrix(NA,nvar,nstart)
  dnorm  <-  c()
  X  <-  list()
  G  <-  list()
  w  <-  list()
  mess <- list()
  for(run in 1:nstart){
    f0 <- fn(x0[,run])
    g0 <- gr(x0[,run])
    
    if((is.na(f0)) || (f0  ==  Inf) || (maxit  == 0) && (prtlevel>0)){
      warning('Gradsamp: function is null or infinite or maxit is zero at initial point')
      f[run]  <-  f0
      x[,run]  <-  x0[,run]
      g[,run]  <-  g0
      dnorm[run]  <-  sqrt(sum(g0*g0))
      X[[run]]  <-  x[,run]
      G[[run]]  <-  g0
      w[[run]]  <-  1
      
    } 
    else {
	tmp  <-  gradsamp1run(fn,gr,x0[,run],nvar,f0,g0,samprad,maxit,normtol,ngrad,fvalquit,prtlevel)
	x[,run]  <-  tmp$x
	f[run]  <-  tmp$f
	g[,run]  <-  tmp$g
	dnorm[run]  <-  tmp$dnorm
	X[[run]]  <-  as.matrix(tmp$X)
	G[[run]]  <-  as.matrix(tmp$G)
	w[[run]]  <-  as.matrix(tmp$w)
	mess[[run]] <- tmp$message
    }
  }
return(list(x  =  x,f  =  f,g  =  g,dnorm  =  dnorm, X  =  X, G  =   G, w  =  w,
	    message = mess))
}
