hanso <-
function(fn,gr=NULL,x0 = NULL,nvar=0,nstart=10,maxit = 1000, normtol = 1e-6, 
		      fvalquit = -Inf, xnormquit = Inf, nvec = 0, prtlevel = 1,
		      strongwolfe = 0, wolfe1 = 1e-4, wolfe2 = 0.5, quitLSfail = 1,
		      ngrad  =  min(100,2*nvar,nvar+10), evaldist = 1e-4, H0 = diag(nvar), scale = 1,
		      samprad  =  c(1e-4,1e-5,1e-6)){
  
  if(!is.null(x0)){
    
    if(class(x0) == "numeric"){
      x0 <- matrix(x0)
      nstart <- 1
      nvar = length(x0)
    }
    else if(class(x0) == "matrix"){
      nvar <- nrow(x0)
      nstart <- ncol(x0)
    }
    else stop("unknown initial value matrix, please enter a numeric vector or matrix")
  }
  else{
    nstart <- 10
    x0 <- matrix(rnorm(nvar*nstart),nvar,nstart)
  }
    
  if(is.null(gr)){
    gr <- function(x){
      grad_nso(fn,x,dir="forward")
    }
  }
  
  tmp <- bfgs(fn,gr,nvar,nstart,x0,maxit, normtol, 
		      fvalquit, xnormquit, nvec , prtlevel ,
		      strongwolfe , wolfe1, wolfe2, quitLSfail,
		      ngrad, evaldist, H0, scale)
  x <- tmp$x
  f <- tmp$f
  d <- tmp$d
  H <- tmp$H
  iter <- tmp$iter
  info <- tmp$info
  X <- tmp$X
  G <- tmp$G
  w <- tmp$w
  
  
  if(length(f)>1){
    indx <- which.min(f)
    f <- f[indx]
    x <- x[,indx]
    d <- d[[indx]]
    H <- H[[indx]]
    X <- X[[indx]]
    G <- G[[indx]]
    w <- w[[indx]]
  }
  else{
    d <- d[[1]]
    H <- H[[1]]
    X <- X[[1]]
    G <- G[[1]]
    w <- w[[1]]
  }
  dnorm <-sqrt(sum(d*d))
  tmp <-postprocess(x,NA,dnorm,X,G,w)
  loc <-tmp$loc
  X <-tmp$X
  G <-tmp$G
  w <-tmp$w
  
  if(is.nainf(f)){
    if(prtlevel>0) warning("Hanso: f is infinite or nan at starting position")
    return(list(x=x,f=f,loc=loc,X=X,G=G,w=w,H=H))
  }
  
  if(f < fvalquit){
    return(list(x=x,f=f,loc=loc,X=X,G=G,w=w,H=H))
  }
  
  if(dnorm < normtol){
    return(list(x=x,f=f,loc=loc,X=X,G=G,w=w,H=H))
  }
  #conditions check
  if(prtlevel) cat("Hanso: Best value found by BFGS = ",f,"\n")
  if(length(samprad)){
    #gradient sampling
    f_BFGS <-f
    dnorm_BFGS <-dnorm
    loc_BFGS <-loc
    d_BFGS <-d
    X_BFGS <-X
    G_BFGS <- G
    w_BFGS <-w
    x0 <-x
    maxit=min(100,maxit)
    nstart <-1
    tmp <-gradsamp(fn,gr,nvar,x0,f0  =  fn(x0), g0  =  gr(x0), samprad,maxit,normtol, ngrad, fvalquit, prtlevel)
    x <-tmp$x
    f <-tmp$f
    g <-tmp$g
    dnorm <-tmp$dnorm
    X <-tmp$X[[1]]
    G <-tmp$G[[1]]
    w <-tmp$w[[1]]
    if(f == f_BFGS){
      if(prtlevel) warning("gradsamp did not reduce f below f_BFGS")
      if(dnorm > dnorm_BFGS){
	loc <-loc_BFGS
	d <-d_BFGS
	X <-X_BFGS
	G <-G_BFGS
	w <-w_BFGS
      }}
      else if(f < f_BFGS) {
	tmp <-postprocess(x,g,dnorm,X,G,w)
	loc <-tmp$loc
	X <-tmp$X
	G <-tmp$G
	w <-tmp$w
      }
      else warning('Hanso: f > f_BFGS, this should not happen')
    }
    cat("Best value found by Gradient Sampling = ",f,"\n")
    return(list(x=x,f=f,loc=loc,X=X,G=G,w=w,H=H))
}
