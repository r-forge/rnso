bfgs1run <- function(fn, gr, nvar, x0, maxit = 1000, normtol = 1e-6, 
		      fvalquit = -Inf, xnormquit = Inf, nvec = 0, prtlevel = 1,
		      strongwolfe = 0, wolfe1 = 1e-4, wolfe2 = 0.5, quitLSfail = 1,
		      ngrad = 2, evaldist = 1e-4, H0 = diag(nvar), scale = 1){
  n <- nvar
  H <- H0
  x <- x0
  f <- fn(x)
  g <- gr(x)
  d <- g
  G <- g
  X <- x
  nG <- 1
  w <- 1
  dnorm <- norm(as.matrix(g),type <- "1")
  if(nvec >0){
    S <- c()
    Y <- c()    
  }
  iter <- 0
  #initializations
  fevalrec <- c()
  xrec <- list() #?
  Hrec <- list() #?
  ##print error msgs
  if(isnaninf(f)) {
      if(prtlevel >0) warning('BFGS: f is infinite or nan at initial iterate')
      info <- 5
      return(list(x = x,f = f,d = d,H = H,iter = iter,info = info,X = X,G = G,
		  w = w,fevalrec = fevalrec,xrec = xrec,Hrec = Hrec))
    }
    
  else if (isnaninf(g)){
      info <- 5
      if (prtlevel>0) warning('BFGS: gradient is infite or nan at initial iterate')
      return(list(x = x,f = f,d = d,H = H,iter = iter,info = info,X = X,G = G,
		  w = w,fevalrec = fevalrec,xrec = xrec,Hrec = Hrec))
  }
  
  else if(dnorm <normtol){
    info <- 0
    if(prtlevel > 0) warning('BFGS: tolerance of gradient satisfied at initial iterate')
    return(list(x = x,f = f,d = d,H = H,iter = iter,info = info,X = X,G = G,
		  w = w,fevalrec = fevalrec,xrec = xrec,Hrec = Hrec))
  }
  
  else if(f < fvalquit){
    info <- 2
    if(prtlevel > 0) warning('BFGS: below target objective at initial iterate')
    return(list(x = x,f = f,d = d,H = H,iter = iter,info = info,X = X,G = G,
		  w = w,fevalrec = fevalrec,xrec = xrec,Hrec = Hrec))
  }
  
  else if(norm(as.matrix(x)) > xnormquit){
    info <- 3
    if(prtlevel > 0) warning('BFGS: norm(x) exceeds limit at initial iterate')
    return(list(x = x,f = f,d = d,H = H,iter = iter,info = info,X = X,G = G,
		  w = w,fevalrec = fevalrec,xrec = xrec,Hrec = Hrec))
  }
  
  
  for(iter in 1:maxit){
    if(nvec == 0){
      p <- -H%*%g
    }
    
    else{
       p <- -1*hgprod(H,g,S,Y) 
    } 
    
    gtp <- t(g)%*%p
    if((gtp >= 0) | is.nan(gtp)){
      if(prtlevel >0) warning("bfgs: non descent direction, quit")
      info <- 6
      return(list(x = x,f = f,d = d,H = H,iter = iter,info = info,X = X,G = G,
		  w = w,fevalrec = fevalrec,xrec = xrec,Hrec = Hrec))
    }
    
    gprev <- g
    if(strongwolfe){ # we do not have linesch_sw()
#       fevalrecline <- nan;
#       tmp <- linesch_sw(x,f,g,p,pars,wolfe1,wolfe2,fvalquit,prtlevel)
#       
#       if(wolfe2 == 0){
#         increase <- 1e-8*(1+alpha)
#         x <- x+increase*p
#       }
#       if(prtlevel>1) warning("exact line sch simulation ")
#       f <- fn(x)
#       g <- gr(x)
    }
    else {
      tmp <- linesch_ww(fn, gr, x, d=p, fn0 = fn(x), gr0 = gr(x),
                        c1 = wolfe1, c2 = wolfe2)
      alpha <- tmp$alpha
      x <- tmp$xalpha
      f <- tmp$falpha
      g <- tmp$galpha
      fail <- tmp$fail
      fevalrecline <- tmp$fevalrec
    }
    
    if(alpha*norm(as.matrix(p))>evaldist){
      nG <- 1
      G <- g
      X <- x
    }
    
    else if(nG<ngrad){
      nG <- nG+1
      G <- cbind(g,G)
      X <- cbind(x,X)
    }
    
    else{
      G <- cbind(g,G[,1:ngrad-1])
      X <- cbind(x,X[,1:ngrad-1])
    }
    
    
    if(nG>1) {
      qtemp <- qpspecial(t(G))
      w <- qtemp$x
      d <- qtemp$d
    }else{
      w <- 1
      d <- g
    }
    
    dnorm <- norm(as.matrix(d))
    
    xrec[[iter]] <- x #xrec is a list
    fevalrec[[iter]] <- fevalrecline
    Hrec[[iter]] <- H
    
#     if(prtlevel>1){  ## no prints needed
#       nfeval <- length(fevalrecline)
#       
#     }
    
    if(f<fvalquit){
      info <- 2
      return(list(x = x,f = f,d = d,H = H,iter = iter,info = info,X = X,G = G,
		  w = w,fevalrec = fevalrec,xrec = xrec,Hrec = Hrec)) 
    }
    
    else if(norm(as.matrix(x))>xnormquit){
      if(prtlevel>0) warning("bfgs: norm exceeds specified limit")
      info <- 3
      return(list(x = x,f = f,d = d,H = H,iter = iter,info = info,X = X,G = G,
		  w = w,fevalrec = fevalrec,xrec = xrec,Hrec = Hrec))
    }
    
    if(fail == 1){
      if(!quitLSfail){
        if(prtlevel>1) warning("bfgs: continue alhough line search failed")
      }else{
        if(prtlevel>0) warning("bfgs: quit at iteration...")
        info <- 7
        return(list(x = x,f = f,d = d,H = H,iter = iter,info = info,X = X,G = G,
		  w = w,fevalrec = fevalrec,xrec = xrec,Hrec = Hrec))
      }
    }
    
    if(fail == -1){
      if(prtlevel>0) warning("bfgs: f may be unbounded below...")
      info <- 8
      return(list(x = x,f = f,d = d,H = H,iter = iter,info = info,X = X,G = G,
		  w <- w,fevalrec = fevalrec,xrec = xrec,Hrec = Hrec))
    }
    if(dnorm <= normtol){
      if(prtlevel>0){
        if(nG == 1) warning("bfgs: gradient norm below tolerance, quite iteration")
        else warning("bfgs: norm of smallest vector in convex hull of gradients below tolerance, quit...")
      }
      info <- 0
      return(list(x = x,f = f,d = d,H = H,iter = iter,info = info,X = X,G = G,
		  w = w,fevalrec = fevalrec,xrec = xrec,Hrec = Hrec))
    }
    
    # ommitted cpu time stuff
    
    s <- alpha*p
    y <- g-gprev
    sty <- t(s)%*%y
    
    if(nvec == 0){
      if(sty > 0){
        if(iter == 1 & scale) H <- as.numeric(sty/(t(y)%*%y))*H
        rho <- as.numeric(1/sty)
        rhoHyst <- rho*(H%*%as.matrix(y))%*%t(s)
        
        H <- H-t(rhoHyst)-rhoHyst+rho*s%*%(t(y)%*%rhoHyst)+rho*s%*%t(s)
      }else{
        if(prtlevel>1) warning("bfgs: sty< <- 0, skipping bfgs update at iteration...")
      }
    }else{
    
      s <- alpha*p
      y <- g-gprev
      if(iter <= nvec){
        S <- cbind(S,s)
        Y <- cbind(Y,y)
      }
      else {
        S <- cbind(S[,2:nvec],s)
        Y <- cbind(Y[,2:nvec],y)
      }
      if(scale) H <- ((t(s)%*%y)/(t(y)%*%y))%*%H0
    }
    
  }
  info <- 1
  return(list(x = x,f = f,d = d,H = H,iter = iter,info = info,X = X,G = G,w = w,fevalrec = fevalrec,xrec = xrec,Hrec = Hrec))
}



















