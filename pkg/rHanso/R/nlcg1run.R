nlcg1run <- function(fn, gr, x0,dir, H0 = NULL, maxit = 1000,  fvalquit = -Inf, 
		      strongwolfe = 1, version = 'C', prtlevel = 0,
                     normtol = 1e-6, xnormquit = Inf, evaldist = 1e-4, ngrad = 0,
                     scale = 1, wolfe1 = 1e-4, wolfe2 = 0.5, quitLSfail = TRUE)
           
           {
		
		is.nainf <- function(x) any(is.na(x) || is.infinite(x))
		frec <- c()
		alpharec <- c()
		x <- x0
		f <- fn(x)
		g <- gr(x)
		iter <- 0
		gnorm <- sqrt(sum(g*g))
		
		if (is.nainf(f) || is.nainf(g)) {
		  
		  mess <- paste("Function or its gradient not defined at initial iterate. iter=",iter)
		  return( list(x = c(x), f = f, g = c(g), frec = frec, alpharec = alpharec,
		  message = mess))
		  }
		if (f < fvalquit) {
		 
		  mess <- paste("Function value found below target at initial iterate.iter = ",iter)
		  return( list(x = c(x), f = f, g = c(g), frec = frec, alpharec = alpharec,
		  message = mess))
		  }
		if (gnorm < normtol) {
		 
		  mess <- paste("Tolerance o gradient satisfied at initial iterate.iter=",iter)
		  return( list(x = c(x), f = f, g = c(g), frec = frec, alpharec = alpharec,
		  message = mess))
		  }  
		
		p <- -g
    
		for(iter in 1:maxit){
		  
		  
		  gtp <- c(t(g) %*% p)
		   if (gtp >= 0 || is.na(gtp)) {
		  errno <- 6
		  mess <- paste("Found non-descent direction only, will quit. iter=",iter)
		  return( list(x = c(x), f = f, g = c(g), frec = frec, alpharec = alpharec,
		  message = mess))
                  
		    }
		  gprev <- g  
		  
		  if(strongwolfe){
		    sls <- linesch_sw(fn, gr, x, dir,d = p, f0 = fn(x), grad0 = gr(x),
                          c1 = wolfe1, c2 = wolfe2, fvalquit, prtlevel)
        alpha <- sls$alpha
		    x <- sls$x
		    f <- sls$f
		    g <- sls$grd
		    fail <- sls$fail
		  }
		  else{		 
		    wls <- linesch_ww(fn, gr, x,dir, d = p, fn0 = fn(x), gr0 = gr(x),
                          c1 = wolfe1, c2 = wolfe2)
		    alpha <- wls$alpha
		    x <- wls$xalpha
		    f <- wls$falpha
		    g <- wls$galpha
		    fail <- wls$fail
		  }
		  
		  gnorm <- sqrt(sum(g*g))
		  frec[iter] <- f
		  alpharec[iter] <- alpha
		  if (f < fvalquit) {
		    errno <- 2
		    mess <- paste("Reached target objective, quit during iteration.",iter)
		    return( list(x = c(x), f = f, g = c(g), frec = frec, alpharec = alpharec,
		  message = mess))
		  }
		  
		  if (fail == 1) {
		    if (!quitLSfail) {
			warning("BFGS: Search continued although line search failed.")
			} else {
			errno <- 7
			mess <- paste("Quit at as line search failed during iteration.",iter)
			return( list(x = c(x), f = f, g = c(g), frec = frec, alpharec = alpharec,
		  message = mess))
			}
			# Function apparently unbounded from below
			} else if (fail == -1) {
			errno <- 8
			mess <- paste("Quit as f may be unbounded from below.",iter)
			return( list(x = c(x), f = f, g = c(g), frec = frec, alpharec = alpharec,
		  message = mess))
			}
			
		if (gnorm <= normtol) {
		  mess <- paste("Gradient norm below tolerance, quit iteration, at iter=",iter)
		  errno <- 0
		  return( list(x = c(x), f = f, g = c(g), frec = frec, alpharec = alpharec,
		  message = mess))
		  }
		
		#update beta according to the version
		
		y <- g-gprev
		if(version == 'P'){
		  nmgprevsq <- sum(gprev*gprev)
		  beta <- (sum(g*y)/nmgprevsq)
		}
		else if(version == 'F'){
		  nmgprevsq <- sum(gprev*gprev)
		  beta <- (sum(g*g)/nmgprevsq)
		}
		else if(version == 'C'){
		  nmgprevsq <- sum(gprev*gprev)
		  beta_pr <- (sum(g*y)/nmgprevsq)
		  beta_fr <- (sum(g*g)/nmgprevsq)
		  if(beta_pr < -beta_fr) beta <- -beta_fr
		  else if(beta_pr > beta_fr) beta <- beta_fr
		  else beta <- beta_pr
		}
		else if(version == 'S'){
		  beta <- (sum(g*y)/sum(p*y))
		}
		else if(version == 'Y'){
		  beta <- (sum(g*g)/sum(p*y))
		}
		else if(version == 'Z'){
		  pty <- sum(p*y)
		  theta <- 2*(sum(y*y)/pty)
		  beta_hz <- (sum((y-theta*p)*g))/(pty)
		  eta <- -1/(sqrt(sum(p*p))*min(0.01,sqrt(sum(gprev*gprev))))
		  beta <- max(beta_hz,eta)
		}
		else if(version == '_'){
		  beta <- 0
		}
		else stop('nlcg: No such version')
		
		
		p <- beta*p - g
		  
		}
		mess='nlcg: number of iterations reached'      
		return( list(x = c(x), f = f, g = c(g), frec = frec, alpharec = alpharec,
			      message = mess))      
		      
		      
           }
